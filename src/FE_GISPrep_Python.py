'''
FastEddy_GISPreprocessor_Python.py
'''

import json
import logging
import math
from pathlib import Path
import re
import sys
import tomllib
import zipfile

import duckdb
import geopandas as gpd
import numpy as np
from pyproj import CRS, Transformer
import rasterio
from rasterio import features
from rasterio.fill import fillnodata
from rasterio.merge import merge
from rasterio.transform import from_bounds as transform_from_bounds
from rasterio.vrt import WarpedVRT
from rasterio.warp import reproject, Resampling, transform_bounds
from rasterio.windows import from_bounds
import requests
from shapely import wkb
import xarray as xr


class FE_GISPrep:

  def __init__(self, setup_params):
    '''Parses parameters and creates projection.'''

    # Parse for shorthand
    self.dname = setup_params['domain_name']
    self.lon_0 = setup_params['lon_0']
    self.lat_0 = setup_params['lat_0']
    self.dh = setup_params['domain_height']
    self.dw = setup_params['domain_width']
    self.cs = setup_params['cell_size']

    # Parse domain bounds (x_min, y_min, x_max, y_max)
    self.db = (
      int(self.dw/-2), int(self.dh/-2), int(self.dw/2), int(self.dh/2)
    )

    # - Specify input url and build output file paths - #

    # Define base_path in the object, too
    self.base_path = Path(setup_params['base_path'])
    
    # Elevation (mandatory)
    self.elev_url = setup_params['elev_url']
    self.elev_tp = self.base_path.joinpath('elev_tiles')
    self.elev_mp = self.base_path.joinpath(f'{self.dname}_elev_mosaic.tif')
    self.elev_fp = self.base_path.joinpath(f'{self.dname}_elev.tif')
    logging.info(f'Elevation data source: {self.elev_url}')
    logging.info(f'Elevation final path: {self.elev_fp}')

    # NLCD (mandatory)
    self.nlcd_url = setup_params['nlcd_url']
    self.nlcd_bp = self.base_path.joinpath(f'{self.dname}_nlcd_buffered.tif')
    self.nlcd_fp = self.base_path.joinpath(f'{self.dname}_nlcd.tif')
    logging.info(f'NLCD data source: {self.nlcd_url}')
    logging.info(f'NLCD final path: {self.nlcd_fp}')

    # Buildings (optional)
    if setup_params['buildings_fp_url']:
      # Building footprints
      self.bldgs_fp_url = setup_params['buildings_fp_url']
      self.bldgs_fp_bp = self.base_path.joinpath(f'{self.dname}_bldgs_fp_buffered.gpkg')
      self.bldgs_fp_fp = self.base_path.joinpath(f'{self.dname}_bldgs_fp_fp.tif')
      logging.info(f'Building footprint data source: {self.bldgs_fp_url}')
      logging.info(f'Building footprints final path: {self.bldgs_fp_fp}')
      
      # LiDAR Point Cloud
      self.lidar_url = setup_params['lidar_url']
      self.lidar_tp = self.base_path.joinpath('lidar_tiles')
      logging.info(f'LiDAR data source: {self.lidar_url}')
      logging.info(f'LiDAR tiles path: {self.lidar_tp}')
    
    else:
      self.bldgs_fp_url = None
      logging.info('No building file data provided.')

    # GIS path is mandatory
    self.nc_fp = self.base_path.joinpath(f'{setup_params['domain_name']}_gis.nc')
    logging.info(f'NetCDF output file: {self.nc_fp}')

    # - #

    # parse proj-string
    self.proj_string = (
      f'+proj=lcc +lon_0={self.lon_0} +lat_0={self.lat_0} '
      f'+lat_1={self.lat_0} +lat_2={self.lat_0}'
    )
    logging.info(f'proj-string: {self.proj_string}')

    return


  def bufferedDomain(self, buffer_size=0.05):
    '''
    Calculates a SW and NE coordinate bounding box from a center point and
    dimensions in meters. The buffer exists to protect against clipped 
    corners during reprojection.

    Parameters
    -----
    buffer_size: Size of the buffer relative to domain (default: 0.05)
    '''

    logging.info('Calculating buffered domain...')
    # Apply the buffer (e.g., 5% means multiplying by 1.05)
    total_x = self.dw * (1.0 + buffer_size)
    total_y = self.dh * (1.0 + buffer_size)
    
    half_x = total_x / 2.0
    half_y = total_y / 2.0
    
    # Approximate meters per degree based on Earth's radius
    meters_per_deg_lat = 111320.0
    meters_per_deg_lon = 111320.0 * math.cos(math.radians(self.lat_0))
    
    # Calculate the offset in degrees
    delta_lat = half_y / meters_per_deg_lat
    delta_lon = half_x / meters_per_deg_lon
    
    self.lat_buff_s = round(self.lat_0 - delta_lat, 4)
    self.lat_buff_n = round(self.lat_0 + delta_lat, 4)
    self.lon_buff_w = round(self.lon_0 - delta_lon, 4)
    self.lon_buff_e = round(self.lon_0 + delta_lon, 4)
    
    logging.info(
      f'Buffered domain: (({self.lat_buff_s}, {self.lon_buff_w}), '
      f'({self.lat_buff_n}, {self.lon_buff_e})).'
    )


  def retrieveBuildingFootprints(self):
    '''
    Fetches Overture Maps building footprints for the domain directly from AWS S3
    using DuckDB, and saves them locally as a GeoPackage.
    '''
    if self.bldgs_fp_bp.is_file():
      logging.info(f'Building footprints already exist at {self.bldgs_fp_bp}.\n')
      return

    logging.info('Connecting to Overture Maps via DuckDB...')
    
    # 1. Initialize DuckDB and load the required extensions
    con = duckdb.connect()
    con.execute('INSTALL spatial; INSTALL httpfs; LOAD spatial; LOAD httpfs;')
    
    # Fetch the latest release string dynamically
    logging.info('Querying Overture STAC catalog for the newest release...')
    catalog_url = 'https://stac.overturemaps.org/catalog.json'
    try:
        latest_release = requests.get(catalog_url).json().get('latest')
        logging.info(f'Latest Overture release identified: {latest_release}')
    except Exception as e:
        logging.error(f'Failed to fetch latest release string: {e}')
        return
    # Build the s3 path
    s3_path = f'{self.bldgs_fp_url}/{latest_release}/theme=buildings/type=building/*'
    
    # 2. The SQL Query
    # Notice the bbox logic: we check for intersection, not just strict inclusion,
    # so we don't accidentally chop off buildings that cross your boundary line.
    query = f'''
      SELECT
        id,
        names.primary AS name,
        height,
        ST_AsWKB(geometry) AS geometry
      FROM read_parquet('{s3_path}', hive_partitioning=1)
      WHERE bbox.xmax >= {self.lon_buff_w}
        AND bbox.xmin <= {self.lon_buff_e}
        AND bbox.ymax >= {self.lat_buff_s}
        AND bbox.ymin <= {self.lat_buff_n}
    '''
    
    logging.info('Executing spatial cloud query (this usually takes 10-20 seconds)...')
    
    # 3. Execute the query and fetch the results as a standard Pandas DataFrame
    df = con.execute(query).df()
    
    if df.empty:
      logging.warning('No buildings found within this domain!')
      return
      
    logging.info(f'Found {len(df)} buildings. Converting to spatial format...')
    
    # 4. Convert the Well-Known Binary (WKB) from DuckDB into Shapely geometries
    df['geometry'] = df['geometry'].apply(
      lambda byte_string: wkb.loads(bytes(byte_string))
    )
    
    # 5. Build the GeoDataFrame and declare the native CRS (Lat/Lon)
    gdf = gpd.GeoDataFrame(df, geometry='geometry', crs='EPSG:4326')
    
    logging.info(f'Saving buildings to: {self.bldgs_fp_bp}')
    
    # 6. Save out as a GeoPackage (.gpkg). It handles mixed data types much better than Shapefiles!
    gdf.to_file(self.bldgs_fp_bp, driver='GPKG')
    
    logging.info('Building footprints downloaded successfully!\n')


  def processBuildings(self):
    '''
    Reads vector building footprints, reprojects them, and burns their height
    values into a raster grid perfectly aligned with the master domain.
    '''
    if self.bldgs_fp_fp.is_file():
      logging.info(f'Rasterized buildings already exist at {self.bldgs_fp_fp}.\n')
      return

    xmin, ymin, xmax, ymax = self.db
    
    # 1. Build the exact same math matrix used for DEM and NLCD
    out_width = int((xmax - xmin) / self.cs)
    out_height = int((ymax - ymin) / self.cs)
    
    out_transform = transform_from_bounds(xmin, ymin, xmax, ymax, out_width, out_height)

    logging.info(f'Loading building footprints from {self.bldgs_fp_bp}...')
    gdf = gpd.read_file(self.bldgs_fp_bp)

    if gdf.empty:
      logging.warning('No buildings to rasterize. Creating empty raster.')
      # Create an empty array of zeros
      burned = np.zeros((out_height, out_width), dtype=np.float32)
    else:
      logging.info(f'Reprojecting buildings to {self.proj_string}...')
      # 2. Reproject the vector data to match your master grid
      gdf = gdf.to_crs(self.proj_string)

      # 3. Handle missing heights
      # Overture might occasionally have a footprint but a null height.
      # We fill NaNs with a default 1-story height (e.g., 4.0 meters) 
      # so the building doesn't disappear from the map.
      gdf['height'] = gdf['height'].fillna(4.0)

      logging.info('Burning building heights into pixel grid...')
      
      # 4. Create an iterable of (geometry, value) pairs for the rasterizer
      shapes = ((geom, value) for geom, value in zip(gdf.geometry, gdf['height']))

      # 5. Rasterize!
      burned = features.rasterize(
        shapes=shapes,
        out_shape=(out_height, out_width),
        transform=out_transform,
        fill=0.0,  # Background pixels get a height of 0
        all_touched=False, # Only burns if the pixel center is inside the building
        dtype=np.float32
      )

    # 6. Save the resulting array to a GeoTIFF
    kwargs = {
      'driver': 'GTiff',
      'crs': self.proj_string,
      'transform': out_transform,
      'width': out_width,
      'height': out_height,
      'count': 1,
      'dtype': 'float32',
      'nodata': 0.0,
      'compress': 'zstd',
      'tiled': True,
      'bigtiff': 'yes'
    }

    logging.info(f'Saving rasterized building heights to: {self.bldgs_fp_fp}')
    with rasterio.open(self.bldgs_fp_fp, 'w', **kwargs) as dst:
      dst.write(burned, 1)

    logging.info('Building rasterization complete!\n')


  def listElev1m(self):
    '''
    Lists all the 1m DEM data tiles that at least partially fall within the
    buffered domain. Returns the list.
    '''
    # TNM Access API expects EPSG:4326 (Lon/Lat)
    # Format: min_lon, min_lat, max_lon, max_lat
    bbox = f'{self.lon_buff_w},{self.lat_buff_s},{self.lon_buff_e},{self.lat_buff_n}'
    
    # Set query parameters
    params = {
      'datasets': 'Digital Elevation Model (DEM) 1 meter',
      'prodFormats': 'GeoTIFF',
      'bbox': bbox,
      'max': 1000,   # Pull massive chunks per page to speed things up
      'offset': 0
    }
      
    elev_urls = []
    logging.info(f'Accessing TNM API for bounding box: {bbox}...\n')
      
    while True:
      response = requests.get(self.elev_url, params=params)
      response.raise_for_status()
      
      try:
        data = response.json()
      except requests.exceptions.JSONDecodeError:
        print('\n[ERROR] The server returned invalid JSON. Here is the raw response:')
        print('-' * 40)
        print(response.text)  # This prints the raw text/HTML the server sent back
        print('-' * 40)
        break
      items = data.get('items', [])
      
      if not items:
        break
              
      logging.info(
        f'Scanning page (offset {params['offset']})... '
        f'found {len(items)} items to inspect.')
          
      for item in items:
        url = item.get('downloadURL', '')
              
        # Filter items
        if '/Elevation/1m/' in url and url.lower().endswith(('.tif', '.tiff')):
          elev_urls.append(url)
                  
      params['offset'] += len(items)
          
      # Break if we've hit the end of the server's records
      if params['offset'] >= data.get("total", 0):
        break
              
    # Remove any duplicates just in case
    elev_urls = list(set(elev_urls))

    # Print number of matching tiles
    logging.info(f'Number of matching tiles: {len(elev_urls)}.')
      
    return elev_urls


  def downloadTiles(self, urls, tile_path):
    '''Downloads the tiles directly to the specified folder using chunking.'''
    tile_path.mkdir(parents=False, exist_ok=True)
      
    for url in urls:
      filename = url.split('/')[-1]
      filepath = tile_path.joinpath(filename)
          
      if filepath.is_file():
        logging.info(f'  Skipping {filename} (Already exists)')
        continue
              
      logging.info(f'  Downloading {filename}...')
      with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(filepath, 'wb') as f:
          for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
                      
    logging.info('\nAll downloads complete!')


  def extractYear(self, filepath):
    '''Searches the filepath for a 4-digit year starting with 20.'''
    match = re.search(r'(20\d{2})', str(filepath))
    # If a year is found, return it as an integer. Otherwise, return 0.
    
    return int(match.group(1)) if match else 0


  def stitchElev(self, input_folder, output_filepath):
    
    # Check for existing file
    if output_filepath.is_file():
      logging.info(f'Elevation mosaic already exists at {output_filepath}.\n')
      return
    
    # Grab all .tif files in the folder
    logging.info(f'Scanning {input_folder} for TIFFs...')
    tif_files = list(input_folder.glob('*.tif'))
    
    if not tif_files:
      logging.warning('No TIFF files found!')
      return
    else:
      logging.info(
        f'{len(tif_files)} tif files found:\n{'\n'.join(f'- {x}' for x in tif_files)}'
      )

    # Sort files from NEWEST to OLDEST; rasterio.merge gap fills missing pixels
    # reverse=True ensures 2020 comes before 2013 in the list.
    tif_files.sort(key=self.extractYear, reverse=True)

    
    # -- Merge Logic -- #
    logging.info(
      '\nVirtually reprojecting and stitching... (This may take a minute...)'
    )
    sources = []
    vrt_list = []
    
    try:
      # Open files and wrap each in a virtual reprojector pointed at target CRS
      for tif in tif_files:
        src = rasterio.open(tif)
        sources.append(src)
        
        vrt = WarpedVRT(src, crs=self.proj_string)
        vrt_list.append(vrt)
        
      # Merge the virtually reprojected tiles!
      mosaic, out_trans = merge(vrt_list)
      
      # Grab metadata from the first VRT to base our output on
      out_meta = vrt_list[0].meta.copy()
      
    finally:
      # Clean up memory by closing all VRTs and source datasets
      for vrt in vrt_list:
        vrt.close()
      for src in sources:
        src.close()
    # -- #

    # Update the metadata for the new stitched grid
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "compress": "zstd",
        'bigtiff': 'yes',
        'tiled': True,
        'num_threads': 'all_cpus'
    })

    # Write the final array to disk
    logging.info(f'\nSaving master DEM to: {output_filepath}')
    with rasterio.open(output_filepath, 'w', **out_meta) as dest:
        dest.write(mosaic)
        
    logging.info('Stitching complete!')


  def prepareNlcdSource(self):
    '''
    Checks if the NLCD parameter is a web URL to a zip file. If so, downloads 
    and extracts it. Otherwise, verifies the local path. Updates the internal 
    NLCD path to point directly to the unzipped raster file.
    '''

    # Create directory to store the NLCD file
    nlcd_dir = self.base_path.joinpath('nlcd_conus')
    nlcd_dir.mkdir(parents=False, exist_ok=True)
    
    # Specify nlcd zip and unzipped file paths
    nlcd_zip_path = nlcd_dir.joinpath(Path(self.nlcd_url).name)
    nlcd_tif_path = nlcd_dir.joinpath(f'{Path(self.nlcd_url).stem}.tif')
    
    # Download if we haven't already
    if not nlcd_zip_path.is_file():
      logging.info(f'Downloading NLCD zip from {self.nlcd_url} (This may take a while)...')
      with requests.get(self.nlcd_url, stream=True) as r:
        r.raise_for_status()
        with open(nlcd_zip_path, 'wb') as f:
          for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
      logging.info('NLCD download complete.')
    else:
      logging.info('NLCD zip archive already exists locally.')
      
    # Unzip if we haven't already
    if not nlcd_tif_path.exists():
      logging.info('Extracting NLCD archive...')
      with zipfile.ZipFile(nlcd_zip_path, 'r') as zip_ref:
        zip_ref.extractall(nlcd_dir)
      logging.info('Extraction complete.')
      
    # Find the raster file inside the extracted folder
    # MRLC uses either .img or .tif
    raster_files = list(nlcd_dir.rglob('*.tif'))
    
    if not raster_files:
      logging.error('No .tif files found in the extracted NLCD folder!')
      sys.exit(1)
      
    # Update the class variable to point to the actual file
    self.nlcd_url = raster_files[0]
    logging.info(f'NLCD source successfully set to: {self.nlcd_url}')


  def bufferedNLCD(self):
    '''
    Clips the CONUS NLCD dataset to the specified bounding box using a 
    memory-efficient windowed read, and saves it to the output path.
    '''
    
    # 1. Coordinate Setup - Use rasterio's built-in bounds transformer
    # This safely calculates the max extents of the curved polygon!
    left, bottom, right, top = transform_bounds(
      'EPSG:4326', 'EPSG:5070', 
      self.lon_buff_w, self.lat_buff_s, self.lon_buff_e, self.lat_buff_n
    )

    logging.info(f'Opening NLCD source: {self.nlcd_url}')
    
    with rasterio.open(self.nlcd_url) as src:
      # 2. Calculate the exact pixel window that corresponds to our bounding box
      window = from_bounds(left, bottom, right, top, src.transform)
      window = window.round_lengths().round_offsets()
      
      logging.info('Reading windowed data into memory...')
      # 3. Read ONLY the data inside our window
      clipped_array = src.read(1, window=window)
      window_transform = src.window_transform(window)
      
      # 4. Prepare the metadata for the new clipped file
      out_meta = src.meta.copy()
      out_meta.update({
        'driver': 'GTiff',
        'height': window.height,
        'width': window.width,
        'transform': window_transform,
        'compress': 'zstd'
      })
      
      logging.info(f'Writing clipped NLCD to: {self.nlcd_bp}')
      # 5. Write the clipped array to the destination path
      with rasterio.open(self.nlcd_bp, 'w', **out_meta) as dest:
        dest.write(clipped_array, 1)
        
    logging.info('NLCD clipping complete!')


  def processRaster(self, input_path, output_path, resampling_method):
    """
    Reprojects, resamples to a specific resolution, and clips to a bounding box.
    """

    # Check for existing file
    if output_path.is_file():
      logging.info(f'Final data already exists at {output_path}.\n')
      return
    
    xmin, ymin, xmax, ymax = self.db
    
    # Calculate the exact number of pixels needed for the new grid.
    out_width = int((xmax - xmin) / self.cs)
    out_height = int((ymax - ymin) / self.cs)
    
    # Build the exact affine transformation matrix for the new clipped grid
    out_transform = transform_from_bounds(xmin, ymin, xmax, ymax, out_width, out_height)
    
    logging.info(f"Processing: {input_path}")
    logging.info(f"Target dimensions: {out_width} x {out_height} pixels")
    
    with rasterio.open(input_path) as src:
      kwargs = src.meta.copy()
      kwargs.update({
        'crs': self.proj_string,
        'transform': out_transform,
        'width': out_width,
        'height': out_height,
        'compress': 'zstd',  # Keep file sizes manageable
        'tiled': True,
        'bigtiff': 'yes',
        'num_threads': 'all_cpus'
      })

      with rasterio.open(output_path, 'w', **kwargs) as dst:
        # Simple loop without the progress bar overhead
        for i in range(1, src.count + 1):
          reproject(
            source=rasterio.band(src, i),
            destination=rasterio.band(dst, i),
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=out_transform,
            dst_crs=self.proj_string,
            # Use Resampling.nearest for categorical masks/flags, bilinear for continuous data
            resampling=resampling_method
          )
          
    logging.info(f"Saved to: {output_path}\n")


  def genCoords(self):
    '''Generates X, Y, LON, and LAT coordinate arrays.'''

    # create x and y indices
    logging.info('Creating x and y indices...')
    self.xs = np.arange(self.db[0], self.db[2], self.cs)
    self.ys = np.arange(self.db[1], self.db[3], self.cs)

    xy = np.meshgrid(self.xs, self.ys)
    logging.info('Index creation complete.')

    # transform to geodetic coordinates
    logging.info('Transforming to WGS 1984 geodetic coordinates...')
    crs = CRS.from_proj4(self.proj_string)
    t = Transformer.from_crs(crs, 4326)  # 4326: WGS 1984
    self.lats, self.lons = t.transform(xy[0], xy[1])
    logging.info('Coordinate transform complete.')

    return


  def gisNetCDF(self):
    '''Creates GIS NetCDF without building heights.'''

    logging.info('Building final NetCDF dataset...')

    # Open mandatory datasets
    with rasterio.open(self.elev_fp) as src:
      elev = src.read(1)
    with rasterio.open(self.nlcd_fp) as src:
      nlcd = src.read(1)

    # Fill missing values in elev and nlcd (if any) with nearest neighbor
    # Assumes missing values set to -9999
    elev = np.where(elev == -9999, np.nan, elev)
    elev = fillnodata(elev, mask=~np.isnan(elev))
    nlcd = np.where(nlcd == -9999, np.nan, nlcd)
    nlcd = fillnodata(nlcd, mask=~np.isnan(nlcd))

    # Flip GIS fields along the y-axis (move origin from NW to SW)
    elev = np.flip(elev, axis=0)
    nlcd = np.flip(nlcd, axis=0)

    # Build base data dictionary
    data_vars = {
      'x': (['x'], self.xs),
      'y': (['y'], self.ys),
      'elevation': (['y', 'x'], elev),
      'lat': (['y', 'x'], self.lats),
      'lon': (['y', 'x'], self.lons),
      'LandCover': (['y', 'x'], nlcd),
      'cellsize': self.cs
    }

    # Conditionally add building data
    if self.bldgs_fp_url:
      with rasterio.open(self.bldgs_fp_fp) as src:
        buildings = src.read(1)
      buildings = np.flip(buildings, axis=0)
      data_vars['BuildingHeights'] = (['y', 'x'], buildings)

    # Build and save Dataset
    ds = xr.Dataset(
      data_vars=data_vars,
      attrs=dict(description='NetCDF file created from automated Python pipeline')
    )
    ds.to_netcdf(self.nc_fp)

    logging.info('Dataset creation complete. Good bye.')

    return


if __name__ == '__main__':
  # - Initialize processing object - #
  if len(sys.argv) == 2:

    # Load parameters
    # NOTE: Will likely switch to .json only.
    param_file = sys.argv[1]
    if Path(param_file).suffix == '.toml':
      with open(param_file, 'rb') as pf:
        setup_params = tomllib.load(pf)
    elif Path(param_file).suffix == '.json':
      with open(param_file, 'r') as pf:
        setup_params = json.load(pf)
      
    # Create base path if it does not exist
    domain_name = setup_params['domain_name']
    base_path = Path(setup_params['base_path'])

    base_path.mkdir(parents=True, exist_ok=True)
    
    # Configure the global logger
    log_file_name = f'{domain_name}_log.txt'
    log_file = base_path.joinpath(log_file_name)
    
    logging.basicConfig(
      filename=log_file,
      filemode='w',
      level=logging.INFO,
      format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Initialize processing object
    logging.info(f'Preparing FastEddy GIS for {domain_name}')
    FGP = FE_GISPrep(setup_params)
  else:
    print(
      'Usage: python -m FE_GISPrep_Python {path_to_parameter_file}'
    )
    sys.exit(1)
  # - #

  # - Create a buffered domain - #
  FGP.bufferedDomain()
  # - #

  # - Retrieve data files - #
  # Elevation
  elev_urls = FGP.listElev1m()
  FGP.downloadTiles(elev_urls, FGP.elev_tp)

  # NLCD (if not yet stored locally)
  nlcd_str = str(FGP.nlcd_url).strip()
  if nlcd_str.startswith('http') and nlcd_str.lower().endswith('.zip'):
    FGP.prepareNlcdSource()
    
  # Buildings
  FGP.retrieveBuildingFootprints()
  # - #

  # - Preprocess input data - #
  # Mosaic elevation data
  FGP.stitchElev(FGP.elev_tp, FGP.elev_mp)
  # Clip NLCD to buffered domain
  FGP.bufferedNLCD()
  # - #
  
  # - Reproject, resample, and clip buffered domain data files - #
  FGP.processRaster(FGP.elev_mp, FGP.elev_fp, Resampling.bilinear)
  FGP.processRaster(FGP.nlcd_bp, FGP.nlcd_fp, Resampling.nearest)
  FGP.processBuildings()
  # - #

  # - Build FastEddy GIS dataset - #
  FGP.genCoords()
  FGP.gisNetCDF()
  # - #