'''
FastEddy_GISPreprocessor_Python.py
'''

import concurrent.futures
import json
import logging
import math
import multiprocessing
from pathlib import Path
import re
import sys
import tomllib
import zipfile

import duckdb
import geopandas as gpd
import laspy
import numpy as np
import pandas as pd
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


# Open-access global and European data sources (no account needed)
COPERNICUS_DEM_URL = 'https://copernicus-dem-30m.s3.amazonaws.com'
WORLDCOVER_URL = 'https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map'
CORINE_URL = 'https://image.discomap.eea.europa.eu/arcgis/rest/services/Corine/CLC2018_LAEA/MapServer/0'
SOILGRIDS_URL = 'https://files.isric.org/soilgrids/latest/data'

LANDCOVER_SOURCES = {
  'nlcd': 'USGS Annual NLCD land cover (CONUS)',
  'worldcover': 'ESA WorldCover 10 m 2021 v200 (Zanaga et al. 2022), CC-BY 4.0',
  'corine': 'CORINE Land Cover 2018 (CLC2018, level-3 codes), Copernicus Land Monitoring Service, EEA',
}
SOIL_NAMES = {'polaris': 'POLARIS', 'soilgrids': 'SoilGrids 2.0'}
SOIL_SOURCES = {
  'polaris': 'POLARIS v1.0, Chaney et al. (2019), Water Resour. Res. 55:2916',
  'soilgrids': 'SoilGrids 2.0, Poggio et al. (2021), SOIL 7:217, CC-BY 4.0',
}


def fallback_utm_crs(lon, lat=None):
  '''
  UTM CRS assumed for LAZ files without a CRS record: NAD83 UTM over North America (zones 1-23),
  ETRS89 UTM over Europe (zones 28-38), WGS 84 UTM elsewhere (north or south by the latitude).
  '''
  zone = math.floor((lon + 180) / 6) + 1
  if lat is None or (lat >= 0 and zone <= 23 and lon < -50):
    return CRS.from_epsg(26900 + zone)
  if 34 <= lat <= 72 and 28 <= zone <= 38:
    return CRS.from_epsg(25800 + zone)
  return CRS.from_epsg((32600 if lat >= 0 else 32700) + zone)


def process_laz_worker(file, proj_string, elev_mp, lon_0, db, cs, height_filters, mode='raster', gdf_bldgs=None,
                       lat_0=None):
  '''
  Standalone worker for extracting building heights from LAZ point clouds.
  Handles both vector (polygon) and raster (gridded) LOD generation.

  Parameters
  ----------
  file : pathlib.Path
    Path to the .laz file.
  proj_string : str
    PROJ string for the master domain CRS.
  elev_mp : pathlib.Path
    Path to the master elevation mosaic for height normalization.
  lon_0 : float
    Central longitude to determine the native UTM zone if missing.
  db : tuple
    Domain bounds (xmin, ymin, xmax, ymax).
  cs : int/float
    Cell size (resolution) in meters.
  height_filters : str
    'strict' or 'loose' filtering logic for canopy/noise removal.
  mode : str, optional
    'raster' (LOD-1) or 'polygon' (LOD-0). Defaults to 'raster'.
  gdf_bldgs : gpd.GeoDataFrame, optional
    Building footprints for spatial joins (Required if mode='polygon').

  Returns
  -------
  numpy.ndarray or pandas.Series or None
    If mode='raster': Returns 2D numpy array [row, col, clean_height].
    If mode='polygon': Returns pandas Series of clean heights indexed by building ID.
    Returns None if no valid building points are found.
  '''
  with laspy.open(file) as f:
    las = f.read()

  unique_classes = np.unique(las.classification)
  if 6 in unique_classes:
    mask = las.classification == 6
  elif 1 in unique_classes:
    mask = las.classification == 1
  else:
    return None

  bldg_pts = las.points[mask]
  if len(bldg_pts) == 0:
    return None

  x_raw, y_raw, z_raw = bldg_pts.x, bldg_pts.y, bldg_pts.z

  master_crs = CRS.from_proj4(proj_string)
  native_crs_wkt = las.vlrs.get('WktCoordinateSystemVlr')

  # --- CRITICAL RAM CLEANUP ---
  del las
  del bldg_pts

  if native_crs_wkt:
    native_crs = CRS.from_wkt(native_crs_wkt[0].string)
  else:
    native_crs = fallback_utm_crs(lon_0, lat_0)

  transformer = Transformer.from_crs(native_crs, master_crs, always_xy=True)
  x_proj, y_proj = transformer.transform(x_raw, y_raw)

  minx, maxx = np.min(x_proj), np.max(x_proj)
  miny, maxy = np.min(y_proj), np.max(y_proj)

  with rasterio.open(elev_mp) as dem_src:
    window = from_bounds(minx, miny, maxx, maxy, dem_src.transform)
    window = window.round_lengths().round_offsets()
    dem_chunk = dem_src.read(1, window=window)
    chunk_transform = dem_src.window_transform(window)

  inv_transform = ~chunk_transform
  cols_dem, rows_dem = inv_transform * (x_proj, y_proj)
  cols_dem = np.clip(np.floor(cols_dem).astype(int), 0, dem_chunk.shape[1] - 1)
  rows_dem = np.clip(np.floor(rows_dem).astype(int), 0, dem_chunk.shape[0] - 1)

  true_heights = z_raw - dem_chunk[rows_dem, cols_dem]
  valid_mask = true_heights > 2.0
  
  if not np.any(valid_mask):
    return None

  x_valid = x_proj[valid_mask]
  y_valid = y_proj[valid_mask]
  h_valid = true_heights[valid_mask]

  # ---------------------------------------------------------
  # MODE: POLYGON (LOD-0)
  # ---------------------------------------------------------
  if mode == 'polygon':
    df = gpd.GeoDataFrame(
      {'height': h_valid}, 
      geometry=gpd.points_from_xy(x_valid, y_valid), 
      crs=proj_string
    )

    joined = gpd.sjoin(df, gdf_bldgs, how='inner', predicate='within')
    if joined.empty:
      return None

    res_df = joined.groupby('id')['height_left'].quantile(0.9).reset_index(name='p90')

    if height_filters == 'strict':
      p75_df = joined.groupby('id')['height_left'].quantile(0.75).reset_index(name='p75')
      res_df['p75'] = p75_df['p75']
      spike_mask = res_df['p90'] > (3.0 * res_df['p75'])
      res_df['clean_height'] = np.where(spike_mask, res_df['p75'], res_df['p90'])
    else:
      res_df['clean_height'] = res_df['p90']

    area_df = joined.groupby('id')['shape_area'].first().reset_index(name='area')
    res_df['area'] = area_df['area']
    small_mask = (res_df['clean_height'] > 24.0) & (res_df['area'] < 500.0)
    res_df['clean_height'] = np.where(small_mask, 8.4, res_df['clean_height'])

    return res_df.set_index('id')['clean_height']

  # ---------------------------------------------------------
  # MODE: RASTER (LOD-1)
  # ---------------------------------------------------------
  elif mode == 'raster':
    xmin, ymin, xmax, ymax = db
    fe_cols = np.floor((x_valid - xmin) / cs).astype(int)
    fe_rows = np.floor((ymax - y_valid) / cs).astype(int)

    max_col, max_row = int((xmax - xmin) / cs), int((ymax - ymin) / cs)
    in_bounds = (fe_cols >= 0) & (fe_cols < max_col) & (fe_rows >= 0) & (fe_rows < max_row)

    if not np.any(in_bounds):
      return None

    df = pd.DataFrame({
      'row': fe_rows[in_bounds], 
      'col': fe_cols[in_bounds], 
      'height': h_valid[in_bounds]
    })
    
    res_df = df.groupby(['row', 'col'])['height'].quantile(0.9).reset_index(name='p90')
    
    if height_filters == 'strict':
      p50_df = df.groupby(['row', 'col'])['height'].quantile(0.5).reset_index(name='p50')
      res_df['p50'] = p50_df['p50']
      spike_mask = res_df['p90'] > (3.0 * res_df['p50'])
      res_df['clean_height'] = np.where(spike_mask, res_df['p50'], res_df['p90'])
    else:
      res_df['clean_height'] = res_df['p90']
    
    return res_df[['row', 'col', 'clean_height']].values


class FE_GISPrep:

  def __init__(self, setup_params):
    '''
    Initializes the FastEddy GIS Preprocessor, parses configuration parameters, 
    and constructs standardized file paths.

    Parameters
    ----------
    setup_params : dict
      Dictionary of configuration settings parsed from the .toml/.json file.
    '''
    self.dname = setup_params['domain_name']
    self.lon_0 = setup_params['lon_0']
    self.lat_0 = setup_params['lat_0']
    self.dh = setup_params['domain_height']
    self.dw = setup_params['domain_width']
    self.cs = setup_params['cell_size']
    self.lod = setup_params.get('lod')
    self.height_filters = setup_params.get('height_filters', 'loose').lower()
    
    logging.info(f'Building Level of Detail (LOD) set to: {self.lod}')
    logging.info(f'Height filter parameters set to: {self.height_filters}')

    # Optional projection definition decoupled from the domain center (e.g. to match
    # a WRF parent: proj_lat_1/proj_lat_2 = TRUELAT1/2, proj_lon_0 = STAND_LON,
    # proj_R = 6370000). Without proj_lon_0 the default LCC tangent at the domain
    # center is used.
    self.custom_proj = 'proj_lon_0' in setup_params
    if self.custom_proj:
      self.proj_lon_0 = setup_params['proj_lon_0']
      self.proj_lat_1 = setup_params.get('proj_lat_1', self.lat_0)
      self.proj_lat_2 = setup_params.get('proj_lat_2', self.proj_lat_1)
      self.proj_lat_0 = setup_params.get('proj_lat_0', self.proj_lat_1)
      self.proj_R = setup_params.get('proj_R')
      self.proj_string = (
        f'+proj=lcc +lon_0={self.proj_lon_0} +lat_0={self.proj_lat_0} '
        f'+lat_1={self.proj_lat_1} +lat_2={self.proj_lat_2} +units=m +no_defs'
      )
      if self.proj_R:
        self.proj_string += f' +R={self.proj_R}'
    else:
      self.proj_lon_0 = self.lon_0
      self.proj_lat_0 = self.proj_lat_1 = self.proj_lat_2 = self.lat_0
      self.proj_R = None
      self.proj_string = (
        f'+proj=lcc +lon_0={self.lon_0} +lat_0={self.lat_0} '
        f'+lat_1={self.lat_0} +lat_2={self.lat_0}'
      )
    logging.info(f'proj-string: {self.proj_string}')

    if self.custom_proj:
      # Center the domain on (lat_0, lon_0) within the projection plane
      t = Transformer.from_crs('EPSG:4326', CRS.from_proj4(self.proj_string), always_xy=True)
      xc, yc = t.transform(self.lon_0, self.lat_0)
      logging.info(f'Domain center in projection coordinates: ({xc:.2f}, {yc:.2f}) m')
      self.db = (
        xc - self.dw/2, yc - self.dh/2, xc + self.dw/2, yc + self.dh/2
      )
    else:
      self.db = (
        int(self.dw/-2), int(self.dh/-2), int(self.dw/2), int(self.dh/2)
      )

    # Elevation product (TNM dataset name and download-URL filter) and resampling
    self.elev_dataset = setup_params.get('elev_dataset', 'Digital Elevation Model (DEM) 1 meter')
    self.elev_url_filter = setup_params.get('elev_url_filter', '/Elevation/1m/')
    self.nlcd_resampling = Resampling[
      setup_params.get('landcover_resampling', setup_params.get('nlcd_resampling', 'nearest'))
    ]
    logging.info(f'Land-cover resampling: {self.nlcd_resampling.name}')

    self.base_path = Path(setup_params['base_path'])

    # Elevation source: 'tnm' (USGS 3DEP via the TNM API, default) or 'copernicus' (Copernicus DEM GLO-30, global,
    # 1 x 1 degree COG tiles on AWS open data)
    self.elev_source = setup_params.get('elev_source', 'tnm').lower()
    if self.elev_source == 'tnm':
      self.elev_url = setup_params['elev_url']
      logging.info(f'Elevation dataset: {self.elev_dataset} (URL filter: {self.elev_url_filter})')
    elif self.elev_source == 'copernicus':
      self.elev_url = setup_params.get('elev_url', COPERNICUS_DEM_URL).rstrip('/')
    else:
      logging.error(f"elev_source = '{self.elev_source}' (use 'tnm' or 'copernicus')")
      sys.exit(1)
    self.elev_tp = self.base_path.joinpath('elev_tiles')
    self.elev_mp = self.base_path.joinpath(f'{self.dname}_elev_mosaic.tif')
    self.elev_fp = self.base_path.joinpath(f'{self.dname}_elev.tif')
    logging.info(f'Elevation data source ({self.elev_source}): {self.elev_url}')
    logging.info(f'Elevation final path: {self.elev_fp}')

    # Land-cover source: 'nlcd' (CONUS, default), 'worldcover' (ESA WorldCover 10 m, global) or 'corine'
    # (CORINE Land Cover 2018, Europe, polygons from the EEA map service). The nlcd_* attributes hold the
    # buffered and final land-cover rasters of whichever source is used.
    self.landcover_source = setup_params.get('landcover_source', 'nlcd').lower()
    if self.landcover_source == 'nlcd':
      self.nlcd_url = setup_params['nlcd_url']
    elif self.landcover_source == 'worldcover':
      self.nlcd_url = setup_params.get('worldcover_url', WORLDCOVER_URL).rstrip('/')
    elif self.landcover_source == 'corine':
      self.nlcd_url = setup_params.get('corine_url', CORINE_URL).rstrip('/')
      self.corine_gp = self.base_path.joinpath(f'{self.dname}_corine_buffered.gpkg')
    else:
      logging.error(f"landcover_source = '{self.landcover_source}' (use 'nlcd', 'worldcover' or 'corine')")
      sys.exit(1)
    self.nlcd_bp = self.base_path.joinpath(f'{self.dname}_{self.landcover_source}_buffered.tif')
    self.nlcd_fp = self.base_path.joinpath(f'{self.dname}_{self.landcover_source}.tif')
    logging.info(f'Land-cover data source ({self.landcover_source}): {self.nlcd_url}')
    logging.info(f'Land-cover final path: {self.nlcd_fp}')

    self.bldgs_fp_url = setup_params.get('buildings_fp_url', '')
    if self.bldgs_fp_url:
      self.bldgs_fp_bp = self.base_path.joinpath(f'{self.dname}_bldgs_fp_buffered.gpkg')
      self.bldgs_fp_fp = self.base_path.joinpath(f'{self.dname}_bldgs_fp_fp.tif')
      logging.info(f'Building footprint data source: {self.bldgs_fp_url}')
      logging.info(f'Building footprints final path: {self.bldgs_fp_fp}')
      
      self.lidar_url = setup_params.get('lidar_url', '')
      self.lidar_tp = self.base_path.joinpath('lidar_tiles')
      logging.info(f'LiDAR data source: {self.lidar_url}')
      logging.info(f'LiDAR tiles path: {self.lidar_tp}')
    else:
      self.bldgs_fp_url = None
      self.lidar_url = None
      logging.info('No building file data provided.')

    # Soil properties per depth interval; optional. soil_source 'polaris' (Chaney et al. 2019, WRR; 30 m, CONUS
    # only; also selected by a non-empty polaris_url) or 'soilgrids' (ISRIC SoilGrids 2.0, Poggio et al. 2021;
    # 250 m, global). Both give sand/silt/clay at the same six depth intervals.
    self.polaris_url = setup_params.get('polaris_url', '')
    self.soil_source = setup_params.get('soil_source', 'polaris' if self.polaris_url else '').lower()
    if self.soil_source:
      if self.soil_source == 'polaris':
        if not self.polaris_url:
          logging.error("soil_source = 'polaris' needs polaris_url")
          sys.exit(1)
        self.polaris_url = self.polaris_url.rstrip('/')
        self.polaris_tp = self.base_path.joinpath('polaris_tiles')
      elif self.soil_source == 'soilgrids':
        self.polaris_url = setup_params.get('soilgrids_url', SOILGRIDS_URL).rstrip('/')
      else:
        logging.error(f"soil_source = '{self.soil_source}' (use 'polaris' or 'soilgrids')")
        sys.exit(1)
      self.polaris_vars = setup_params.get('soil_vars', setup_params.get('polaris_vars', ['sand', 'silt', 'clay']))
      self.polaris_stat = setup_params.get('soil_stat', setup_params.get('polaris_stat', 'mean'))
      self.polaris_depths = setup_params.get('soil_depths', setup_params.get(
        'polaris_depths', ['0_5', '5_15', '15_30', '30_60', '60_100', '100_200']))
      logging.info(f'Soil data source ({self.soil_source}): {self.polaris_url} ({self.polaris_stat} of '
                   f'{self.polaris_vars}, depths {self.polaris_depths} cm)')
    else:
      self.polaris_url = None
      logging.info('No soil data requested.')

    self.nc_fp = self.base_path.joinpath(f'{self.dname}_gis.nc')
    logging.info(f'NetCDF output file: {self.nc_fp}')

    # Initialize empty lists to store targeted file paths
    self.elev_files = []
    self.laz_files = []


  def bufferedDomain(self, buffer_size=0.05):
    '''
    Calculates a SW and NE coordinate bounding box from a center point and
    dimensions in meters. The buffer exists to protect against clipped 
    corners during reprojection.

    Parameters
    ----------
    buffer_size : float, optional
      Size of the spatial buffer relative to the domain (default is 0.05 / 5%).
    '''
    logging.info('Calculating buffered domain...')
    if self.custom_proj:
      # Transform the buffered projected box (densified edges) to geographic bounds;
      # a center/cos(lat) estimate clips corners of large or off-origin domains.
      bx = 0.5 * buffer_size * self.dw
      by = 0.5 * buffer_size * self.dh
      w, s, e, n = transform_bounds(
        self.proj_string, 'EPSG:4326',
        self.db[0] - bx, self.db[1] - by, self.db[2] + bx, self.db[3] + by,
        densify_pts=101
      )
      self.lat_buff_s = math.floor(s * 1e4) / 1e4
      self.lat_buff_n = math.ceil(n * 1e4) / 1e4
      self.lon_buff_w = math.floor(w * 1e4) / 1e4
      self.lon_buff_e = math.ceil(e * 1e4) / 1e4
      logging.info(
        f'Buffered domain: (({self.lat_buff_s}, {self.lon_buff_w}), '
        f'({self.lat_buff_n}, {self.lon_buff_e})).\n'
      )
      return

    total_x = self.dw * (1.0 + buffer_size)
    total_y = self.dh * (1.0 + buffer_size)
    
    half_x = total_x / 2.0
    half_y = total_y / 2.0
    
    meters_per_deg_lat = 111320.0
    meters_per_deg_lon = 111320.0 * math.cos(math.radians(self.lat_0))
    
    delta_lat = half_y / meters_per_deg_lat
    delta_lon = half_x / meters_per_deg_lon
    
    self.lat_buff_s = round(self.lat_0 - delta_lat, 4)
    self.lat_buff_n = round(self.lat_0 + delta_lat, 4)
    self.lon_buff_w = round(self.lon_0 - delta_lon, 4)
    self.lon_buff_e = round(self.lon_0 + delta_lon, 4)
    
    logging.info(
      f'Buffered domain: (({self.lat_buff_s}, {self.lon_buff_w}), '
      f'({self.lat_buff_n}, {self.lon_buff_e})).\n'
    )


  def retrieveBuildingFootprints(self):
    '''
    Fetches Overture Maps building footprints directly from AWS S3 using 
    DuckDB spatial extensions, and saves them locally as a GeoPackage.
    '''
    if self.bldgs_fp_bp.is_file():
      logging.info(f'Building footprints already exist at {self.bldgs_fp_bp}.\n')
      return

    logging.info('Connecting to Overture Maps via DuckDB...')
    
    con = duckdb.connect()
    con.execute('INSTALL spatial; INSTALL httpfs; LOAD spatial; LOAD httpfs;')
    
    catalog_url = 'https://stac.overturemaps.org/catalog.json'
    try:
      latest_release = requests.get(catalog_url).json().get('latest')
      logging.info(f'Latest Overture release identified: {latest_release}')
    except Exception as e:
      logging.error(f'Failed to fetch latest release string: {e}')
      return
      
    s3_path = f'{self.bldgs_fp_url}/{latest_release}/theme=buildings/type=building/*'
    
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
    df = con.execute(query).df()
    
    if df.empty:
      logging.warning('No buildings found within this domain!')
      return
      
    logging.info(f'Found {len(df)} buildings. Converting to spatial format...')
    df['geometry'] = df['geometry'].apply(
      lambda byte_string: wkb.loads(bytes(byte_string))
    )
    
    gdf = gpd.GeoDataFrame(df, geometry='geometry', crs='EPSG:4326')
    
    logging.info(f'Saving buildings to: {self.bldgs_fp_bp}')
    gdf.to_file(self.bldgs_fp_bp, driver='GPKG')
    logging.info('Building footprints downloaded successfully!\n')


  def extractLidarHeights(self):
    '''
    Extracts high-resolution roof heights from LiDAR files in parallel, applying 
    noise filters, and maps the resulting maximum height back to the vector footprints.
    '''
    logging.info('Extracting LiDAR building heights (Memory-Safe Multiprocessing)...')
    
    laz_files = self.laz_files
    if not laz_files:
      logging.warning('No LAZ files found in directory.')
      return
      
    if not self.bldgs_fp_bp.is_file():
      logging.error('Building footprints not found! Cannot extract LiDAR.')
      return

    logging.info('Loading master building footprints into RAM...')
    gdf_bldgs = gpd.read_file(self.bldgs_fp_bp)
    gdf_bldgs_proj = gdf_bldgs.to_crs(self.proj_string)
    gdf_bldgs_proj['shape_area'] = gdf_bldgs_proj.geometry.area

    all_p90_series = []
    max_workers = min(2, max(1, int(multiprocessing.cpu_count() * 0.6)))
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
      future_to_file = {
        executor.submit(
          process_laz_worker, file, self.proj_string, self.elev_mp, self.lon_0, 
          self.db, self.cs, self.height_filters, 'polygon', gdf_bldgs_proj, lat_0=self.lat_0
        ): file for file in laz_files
      }
      
      for future in concurrent.futures.as_completed(future_to_file):
        file = future_to_file[future]
        try:
          p90_heights = future.result()
          if p90_heights is not None and not p90_heights.empty:
            all_p90_series.append(p90_heights)
          logging.info(f'  Finished {file.name}')
        except Exception as exc:
          logging.error(f'  {file.name} generated an exception: {exc}')

    if not all_p90_series:
      logging.warning('No LiDAR points intersected with buildings. Using default heights.')
      return

    logging.info('Updating master building footprints with new heights...')
    master_series = pd.concat(all_p90_series)
    final_lidar_heights = master_series.groupby(level=0).max()

    gdf_bldgs['overture_height'] = gdf_bldgs['height'].fillna(4.0)
    gdf_bldgs['lidar_height'] = gdf_bldgs['id'].map(final_lidar_heights)
    gdf_bldgs['height'] = gdf_bldgs['lidar_height'].fillna(gdf_bldgs['overture_height'])
    gdf_bldgs = gdf_bldgs.drop(columns=['overture_height', 'lidar_height'])

    logging.info('Saving updated building heights back to GeoPackage...')
    gdf_bldgs.to_file(self.bldgs_fp_bp, driver='GPKG')
    logging.info('Building height update complete!\n')


  def processBuildings(self):
    '''
    Reads vector footprints and rasterizes them into the master spatial grid.
    If self.lod == 1, dynamically extracts and overlays high-resolution 
    LiDAR pixels into the grid using multiprocessing.
    '''
    xmin, ymin, xmax, ymax = self.db
    out_width = int((xmax - xmin) / self.cs)
    out_height = int((ymax - ymin) / self.cs)

    target_fp = self.bldgs_fp_fp
    if self.lod == 1:
      self.bldgs_lod1_fp = self.base_path.joinpath(f'{self.dname}_bldgs_lod1_fp.tif')
      target_fp = self.bldgs_lod1_fp

    if target_fp.is_file():
      with rasterio.open(target_fp) as chk:
        if chk.width == out_width and chk.height == out_height:
          logging.info(f'Rasterized LOD-{self.lod} buildings already exist at {target_fp}.\n')
          return
        else:
          logging.info(f'Stale dimensions found in {target_fp}. Overwriting...')

    logging.info(f'Generating LOD-{self.lod} building raster...')
    out_transform = transform_from_bounds(xmin, ymin, xmax, ymax, out_width, out_height)

    logging.info('Rasterizing footprint polygons for base canvas...')
    gdf = gpd.read_file(self.bldgs_fp_bp)
    
    if gdf.empty:
      logging.warning('No buildings to rasterize. Creating empty raster.')
      final_array = np.zeros((out_height, out_width), dtype=np.float32)
    else:
      gdf = gdf.to_crs(self.proj_string)
      gdf['shape_area'] = gdf.geometry.area
      gdf['height'] = gdf['height'].fillna(4.0)
      
      small_mask_base = (gdf['height'] > 24.0) & (gdf['shape_area'] < 400.0)
      gdf.loc[small_mask_base, 'height'] = 8.4
      
      shapes = ((geom, value) for geom, value in zip(gdf.geometry, gdf['height']))
      
      base_canvas = features.rasterize(
        shapes=shapes,
        out_shape=(out_height, out_width),
        transform=out_transform,
        fill=0.0,
        all_touched=False,
        dtype=np.float32
      )
      
      if self.lod in [0, None]:
        final_array = base_canvas
        
      elif self.lod == 1:
        laz_files = self.laz_files
        
        if not laz_files:
          logging.warning('No LAZ files found! Falling back to LOD-0.')
          final_array = base_canvas
        else:
          logging.info('Extracting high-res LOD-1 pixels via multiprocessing...')
          footprint_mask = base_canvas > 0
          master_lod1_array = np.zeros((out_height, out_width), dtype=np.float32)
          max_workers = min(2, max(1, int(multiprocessing.cpu_count() * 0.6)))
          
          with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_file = {
              executor.submit(
                process_laz_worker, file, self.proj_string, self.elev_mp, self.lon_0, 
                self.db, self.cs, self.height_filters, 'raster', lat_0=self.lat_0
              ): file for file in laz_files
            }
            
            for future in concurrent.futures.as_completed(future_to_file):
              try:
                pixel_data = future.result()
                if pixel_data is not None:
                  rows = pixel_data[:, 0].astype(int)
                  cols = pixel_data[:, 1].astype(int)
                  heights = pixel_data[:, 2]
                  master_lod1_array[rows, cols] = np.maximum(master_lod1_array[rows, cols], heights)
              except Exception as exc:
                logging.error(f'  LOD-1 worker generated an exception: {exc}')

          logging.info('Applying mask and filters to missing LOD-1 pixels...')
          clean_lod1 = np.where(footprint_mask, master_lod1_array, 0.0)
          
          if self.height_filters == 'strict':
            lod1_spike_mask = (clean_lod1 > (3.0 * base_canvas - 10.0))
            clean_lod1 = np.where(lod1_spike_mask, base_canvas, clean_lod1)
          
          area_shapes = ((geom, value) for geom, value in zip(gdf.geometry, gdf['shape_area']))
          area_canvas = features.rasterize(
            shapes=area_shapes,
            out_shape=(out_height, out_width),
            transform=out_transform,
            fill=0.0,
            all_touched=False,
            dtype=np.float32
          )
          
          lod1_small_mask = (clean_lod1 > 24.0) & (area_canvas > 0) & (area_canvas < 500.0)
          clean_lod1 = np.where(lod1_small_mask, 8.4, clean_lod1)
          
          missing_data_mask = footprint_mask & (clean_lod1 == 0)
          final_lod1 = np.where(missing_data_mask, base_canvas, clean_lod1)
          
          valid_data_mask = final_lod1 > 0
          filled_lod1 = fillnodata(final_lod1, mask=valid_data_mask)
          
          final_array = np.where(footprint_mask, filled_lod1, 0.0).astype(np.float32)

    kwargs = {
      'driver': 'GTiff', 'crs': self.proj_string, 'transform': out_transform,
      'width': out_width, 'height': out_height, 'count': 1,
      'dtype': 'float32', 'nodata': 0.0, 'compress': 'zstd',
      'tiled': True, 'bigtiff': 'yes'
    }

    logging.info(f'Saving rasterized building heights to: {target_fp}')
    with rasterio.open(target_fp, 'w', **kwargs) as dst:
      dst.write(final_array, 1)

    logging.info('Building rasterization complete!\n')


  def query_tnm_api(self, dataset_name, format_name):
    '''
    Queries the USGS National Map API for available tiles intersecting the domain.

    Parameters
    ----------
    dataset_name : str
      The official TNM string for the dataset (e.g., 'Lidar Point Cloud (LPC)').
    format_name : str
      The product format to request (e.g., 'LAZ' or 'GeoTIFF').

    Returns
    -------
    list
      A deduplicated list of download URLs for the requested tiles.
    '''
    bbox = f'{self.lon_buff_w},{self.lat_buff_s},{self.lon_buff_e},{self.lat_buff_n}'
    
    params = {
      'datasets': dataset_name,
      'prodFormats': format_name,
      'bbox': bbox,
      'max': 1000,
      'offset': 0
    }
      
    urls = []
    # Determine the correct base URL based on requested format
    base_url = self.elev_url if format_name == 'GeoTIFF' else self.lidar_url
    
    logging.info(f'Accessing TNM API for bounding box: {bbox}...\n')
      
    while True:
      response = requests.get(base_url, params=params)
      response.raise_for_status()
      
      try:
        data = response.json()
      except requests.exceptions.JSONDecodeError:
        logging.error('The server returned invalid JSON. Raw response:')
        logging.error(response.text)
        break
        
      items = data.get('items', [])
      if not items:
        break
              
      logging.info(f'Scanning page (offset {params["offset"]})... found {len(items)} items.')
          
      for item in items:
        url = item.get('downloadURL', '')
        # Only add valid formats (LiDAR ends with laz, DEM requires the elev_url_filter, default /1m/)
        if format_name == 'LAZ' and url.lower().endswith('.laz'):
          urls.append(url)
        elif format_name == 'GeoTIFF' and self.elev_url_filter in url and url.lower().endswith(('.tif', '.tiff')):
          urls.append(url)
                  
      params['offset'] += len(items)
      if params['offset'] >= data.get('total', 0):
        break
              
    urls = list(set(urls))
    if format_name == 'GeoTIFF':
      urls = self.newestSeamlessTiles(urls)
    logging.info(f'Number of matching tiles: {len(urls)}.')
      
    return urls


  def newestSeamlessTiles(self, urls):
    '''
    Keeps only the newest version of each seamless 1x1 degree DEM tile
    (e.g. USGS_1_n45w099_20220105.tif; an undated name is the current release).
    Other products (e.g. 1 m project tiles) are returned unchanged.

    Parameters
    ----------
    urls : list
      Download URLs returned by the TNM API.

    Returns
    -------
    list
      URLs with superseded seamless-tile versions removed.
    '''
    pattern = re.compile(r'(USGS_\d+_[ns]\d+[ew]\d+)(?:_(\d{8}))?\.tiff?$', re.IGNORECASE)
    newest = {}
    other = []
    for url in urls:
      match = pattern.search(url.split('/')[-1])
      if not match:
        other.append(url)
        continue
      tile, date = match.group(1), match.group(2) or '99999999'
      if tile not in newest or date > newest[tile][0]:
        newest[tile] = (date, url)

    if newest:
      logging.info(f'Seamless DEM: kept newest of {len(urls) - len(other)} versions for {len(newest)} tiles.')
    return other + [url for _, url in newest.values()]


  def downloadTiles(self, urls, tile_path):
    '''
    Downloads remote files directly to the specified local folder using stream chunking.

    Parameters
    ----------
    urls : list
      List of direct download URLs.
    tile_path : pathlib.Path
      Local directory path to save the files.
    '''
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
    '''
    Extracts a 4-digit year starting with '20' from a string or Path.

    Parameters
    ----------
    filepath : str or pathlib.Path
      The string or path to search.

    Returns
    -------
    int
      The extracted year as an integer, or 0 if no year is found.
    '''
    match = re.search(r'(20\d{2})', str(filepath))
    return int(match.group(1)) if match else 0


  def stitchElev(self, output_filepath):
    '''
    Creates a single virtually-stitched mosaic from multiple DEM tiles, prioritizing
    newer tiles when filling gaps, and saves the output to a compressed GeoTIFF.

    Parameters
    ----------
    output_filepath : pathlib.Path
      The desired path for the final stitched output file.
    '''
    if output_filepath.is_file():
      logging.info(f'Elevation mosaic already exists at {output_filepath}.\n')
      return

    tif_files = self.elev_files
    
    if not tif_files:
      logging.warning('No TIFF files found!')
      return
    else:
      logging.info(f'{len(tif_files)} tif files found.')

    tif_files.sort(key=self.extractYear, reverse=True)

    logging.info('\nVirtually reprojecting and stitching... (This may take a minute...)')
    sources = []
    vrt_list = []
    
    try:
      for tif in tif_files:
        src = rasterio.open(tif)
        sources.append(src)
        vrt = WarpedVRT(src, crs=self.proj_string)
        vrt_list.append(vrt)
        
      mosaic, out_trans = merge(vrt_list)
      out_meta = vrt_list[0].meta.copy()
      
    finally:
      for vrt in vrt_list:
        vrt.close()
      for src in sources:
        src.close()

    out_meta.update({
        'driver': 'GTiff',
        'height': mosaic.shape[1],
        'width': mosaic.shape[2],
        'transform': out_trans,
        'compress': 'zstd',
        'bigtiff': 'yes',
        'tiled': True,
        'num_threads': 'all_cpus'
    })

    logging.info(f'\nSaving master DEM to: {output_filepath}')
    with rasterio.open(output_filepath, 'w', **out_meta) as dest:
        dest.write(mosaic)
        
    logging.info('Stitching complete!')


  def prepareNlcdSource(self):
    '''
    Downloads and extracts a zipped NLCD archive if provided as a URL,
    and updates the internal variable to point directly to the local raster file.
    '''
    nlcd_dir = self.base_path.joinpath('nlcd_conus')
    nlcd_dir.mkdir(parents=False, exist_ok=True)
    
    nlcd_zip_path = nlcd_dir.joinpath(Path(self.nlcd_url).name)
    nlcd_tif_path = nlcd_dir.joinpath(f'{Path(self.nlcd_url).stem}.tif')
    
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
      
    if not nlcd_tif_path.exists():
      logging.info('Extracting NLCD archive...')
      with zipfile.ZipFile(nlcd_zip_path, 'r') as zip_ref:
        zip_ref.extractall(nlcd_dir)
      logging.info('Extraction complete.')
      
    raster_files = list(nlcd_dir.rglob('*.tif'))
    
    if not raster_files:
      logging.error('No .tif files found in the extracted NLCD folder!')
      sys.exit(1)
      
    self.nlcd_url = raster_files[0]
    logging.info(f'NLCD source successfully set to: {self.nlcd_url}')


  def bufferedNLCD(self):
    '''
    Performs a memory-efficient windowed read of the national NLCD dataset,
    extracting only the pixels within the buffered domain bounds.
    '''
    logging.info(f'Opening NLCD source: {self.nlcd_url}')
    
    with rasterio.open(self.nlcd_url) as src:
      left, bottom, right, top = transform_bounds(
        'EPSG:4326', src.crs,
        self.lon_buff_w, self.lat_buff_s, self.lon_buff_e, self.lat_buff_n
      )
      window = from_bounds(left, bottom, right, top, src.transform)
      window = window.round_lengths().round_offsets()
      
      logging.info('Reading windowed data into memory...')
      clipped_array = src.read(1, window=window)
      window_transform = src.window_transform(window)
      
      out_meta = src.meta.copy()
      out_meta.update({
        'driver': 'GTiff',
        'height': window.height,
        'width': window.width,
        'transform': window_transform,
        'compress': 'zstd'
      })
      
      logging.info(f'Writing clipped NLCD to: {self.nlcd_bp}')
      with rasterio.open(self.nlcd_bp, 'w', **out_meta) as dest:
        dest.write(clipped_array, 1)
        
    logging.info('NLCD clipping complete!')


  def remoteExists(self, url):
    '''
    True if a remote file exists (HTTP HEAD 200). Global tile sets have no tiles over open ocean.
    '''
    try:
      return requests.head(url, timeout=60, allow_redirects=True).status_code == 200
    except requests.RequestException:
      return False


  def copernicusTiles(self):
    '''
    Lists the Copernicus DEM GLO-30 1 x 1 degree COG tiles (EPSG:4326) that cover the buffered domain and exist
    (there are none over open ocean).

    Returns
    -------
    list of str
      Tile URLs such as .../Copernicus_DSM_COG_10_N37_00_W002_00_DEM/Copernicus_DSM_COG_10_N37_00_W002_00_DEM.tif.
    '''
    urls = []
    for la in range(math.floor(self.lat_buff_s), math.ceil(self.lat_buff_n)):
      for lo in range(math.floor(self.lon_buff_w), math.ceil(self.lon_buff_e)):
        name = (f"Copernicus_DSM_COG_10_{'N' if la >= 0 else 'S'}{abs(la):02d}_00_"
                f"{'E' if lo >= 0 else 'W'}{abs(lo):03d}_00_DEM")
        url = f'{self.elev_url}/{name}/{name}.tif'
        if self.remoteExists(url):
          urls.append(url)
        else:
          logging.info(f'  No Copernicus DEM tile {name} (open ocean)')
    return urls


  def worldcoverTiles(self):
    '''
    Lists the ESA WorldCover 3 x 3 degree tiles (EPSG:4326, named by their SW corner) that cover the buffered
    domain and exist.
    '''
    urls = []
    for la in range(3 * math.floor(self.lat_buff_s / 3), math.ceil(self.lat_buff_n), 3):
      for lo in range(3 * math.floor(self.lon_buff_w / 3), math.ceil(self.lon_buff_e), 3):
        name = (f"ESA_WorldCover_10m_2021_v200_{'N' if la >= 0 else 'S'}{abs(la):02d}"
                f"{'E' if lo >= 0 else 'W'}{abs(lo):03d}_Map.tif")
        url = f'{self.nlcd_url}/{name}'
        if self.remoteExists(url):
          urls.append(url)
        else:
          logging.info(f'  No WorldCover tile {name} (open ocean)')
    return urls


  def bufferedLandCover(self):
    '''
    Writes the land cover of the buffered domain to nlcd_bp (NLCD, WorldCover), or the CORINE polygons of the
    buffered domain to corine_gp.
    '''
    if self.landcover_source == 'nlcd':
      self.bufferedNLCD()
    elif self.landcover_source == 'worldcover':
      self.bufferedWorldCover()
    elif self.landcover_source == 'corine':
      self.bufferedCorine()


  def bufferedWorldCover(self):
    '''
    Reads only the buffered domain from the WorldCover COG tiles (HTTP range requests, no full-tile download)
    and writes it as one GeoTIFF (EPSG:4326, 10 m).
    '''
    if self.nlcd_bp.is_file():
      logging.info(f'WorldCover buffered domain already exists at {self.nlcd_bp}.')
      return
    urls = self.worldcoverTiles()
    if not urls:
      logging.error('No WorldCover tiles cover the domain.')
      sys.exit(1)
    sources = [rasterio.open(f'/vsicurl/{u}') for u in urls]
    try:
      bounds = (self.lon_buff_w, self.lat_buff_s, self.lon_buff_e, self.lat_buff_n)
      mosaic, out_trans = merge(sources, bounds=bounds)
      out_meta = sources[0].meta.copy()
    finally:
      for src in sources:
        src.close()
    out_meta.update({'driver': 'GTiff', 'height': mosaic.shape[1], 'width': mosaic.shape[2],
                     'transform': out_trans, 'compress': 'zstd', 'tiled': True})
    with rasterio.open(self.nlcd_bp, 'w', **out_meta) as dest:
      dest.write(mosaic)
    logging.info(f'WorldCover buffered domain written: {self.nlcd_bp}')


  def bufferedCorine(self):
    '''
    Queries the CORINE Land Cover 2018 polygons that intersect the buffered domain from the EEA map service
    (GeoJSON, paged) and saves them, with their level-3 codes (Code_18), to corine_gp.
    '''
    if self.corine_gp.is_file():
      logging.info(f'CORINE polygons already exist at {self.corine_gp}.')
      return
    params = {
      'geometry': f'{self.lon_buff_w},{self.lat_buff_s},{self.lon_buff_e},{self.lat_buff_n}',
      'geometryType': 'esriGeometryEnvelope', 'inSR': 4326, 'spatialRel': 'esriSpatialRelIntersects',
      'outFields': 'Code_18', 'returnGeometry': 'true', 'outSR': 4326, 'f': 'geojson',
      'orderByFields': 'OBJECTID', 'resultRecordCount': 1000,
    }
    features_all = []
    offset = 0
    while True:
      r = requests.get(f'{self.nlcd_url}/query', params=dict(params, resultOffset=offset), timeout=300)
      r.raise_for_status()
      page = r.json().get('features', [])
      features_all.extend(page)
      logging.info(f'  CORINE polygons {offset} to {offset + len(page)}')
      if len(page) < params['resultRecordCount']:
        break
      offset += len(page)
    if not features_all:
      logging.error('No CORINE polygons in the domain (CORINE covers the EEA member and cooperating countries).')
      sys.exit(1)
    gdf = gpd.GeoDataFrame.from_features(features_all, crs='EPSG:4326')
    gdf['Code_18'] = gdf['Code_18'].astype(int)
    gdf.to_file(self.corine_gp, driver='GPKG')
    logging.info(f'{len(gdf)} CORINE polygons written: {self.corine_gp}')


  def rasterizeCorine(self):
    '''
    Rasterizes the CORINE polygons onto the FastEddy master grid (the class of the polygon containing each cell
    centre), written to nlcd_fp. Cells outside every polygon get 0 (no data).
    '''
    xmin, ymin, xmax, ymax = self.db
    out_width = int((xmax - xmin) / self.cs)
    out_height = int((ymax - ymin) / self.cs)
    if self.nlcd_fp.is_file():
      with rasterio.open(self.nlcd_fp) as chk:
        if chk.width == out_width and chk.height == out_height:
          logging.info(f'Final data already exists and matches target dimensions at {self.nlcd_fp}.\n')
          return
    out_transform = transform_from_bounds(xmin, ymin, xmax, ymax, out_width, out_height)
    gdf = gpd.read_file(self.corine_gp).to_crs(self.proj_string)
    arr = features.rasterize(zip(gdf.geometry, gdf['Code_18']), out_shape=(out_height, out_width),
                             transform=out_transform, fill=0, dtype='int16')
    with rasterio.open(self.nlcd_fp, 'w', driver='GTiff', width=out_width, height=out_height, count=1,
                       dtype='int16', nodata=0, crs=self.proj_string, transform=out_transform,
                       compress='zstd', tiled=True) as dst:
      dst.write(arr, 1)
    logging.info(f'CORINE rasterized onto the master grid: {self.nlcd_fp}\n')


  def soilgridsMosaic(self, var, depth, mosaic_fp):
    '''
    Reads the buffered domain of one SoilGrids 2.0 variable and depth from the ISRIC virtual raster
    (Homolosine, 250 m, HTTP range requests) and writes it as a GeoTIFF in the SoilGrids CRS.
    Texture fractions are stored in g/kg; no data (e.g. water, built-up areas) is -32768.
    '''
    d = depth.replace('_', '-')
    url = f'/vsicurl/{self.polaris_url}/{var}/{var}_{d}cm_{self.polaris_stat}.vrt'
    with rasterio.open(url) as src:
      left, bottom, right, top = transform_bounds(
        'EPSG:4326', src.crs, self.lon_buff_w, self.lat_buff_s, self.lon_buff_e, self.lat_buff_n, densify_pts=101
      )
      window = from_bounds(left, bottom, right, top, src.transform).round_lengths().round_offsets()
      arr = src.read(1, window=window)
      out_meta = src.meta.copy()
      out_meta.update({'driver': 'GTiff', 'height': arr.shape[0], 'width': arr.shape[1],
                       'transform': src.window_transform(window), 'compress': 'zstd'})
    with rasterio.open(mosaic_fp, 'w', **out_meta) as dest:
      dest.write(arr, 1)
    logging.info(f'SoilGrids {var} {depth} cm written: {mosaic_fp}')


  def polarisTiles(self):
    '''
    Lists the POLARIS 1 x 1 degree tiles (EPSG:4326) that cover the buffered domain.

    Returns
    -------
    list of str
      Tile names such as 'lat4041_lon-106-105'.
    '''
    lats = range(math.floor(self.lat_buff_s), math.ceil(self.lat_buff_n))
    lons = range(math.floor(self.lon_buff_w), math.ceil(self.lon_buff_e))
    return [f'lat{la}{la + 1}_lon{lo}{lo + 1}' for la in lats for lo in lons]


  def polarisPath(self, var, depth, kind):
    '''
    Local path of a soil product (POLARIS or SoilGrids) for one variable and depth: kind 'mosaic' (source grid)
    or 'final' (on the FastEddy grid).
    '''
    return self.base_path.joinpath(f'{self.dname}_{self.soil_source}_{var}_{self.polaris_stat}_{depth}_{kind}.tif')


  def downloadPolaris(self):
    '''
    Downloads the POLARIS tiles of every requested variable and depth that cover the buffered domain,
    skipping tiles and products that already exist.
    '''
    tiles = self.polarisTiles()
    logging.info(f'POLARIS tiles covering the domain: {tiles}')
    for var in self.polaris_vars:
      for depth in self.polaris_depths:
        if self.polarisPath(var, depth, 'final').is_file() or self.polarisPath(var, depth, 'mosaic').is_file():
          continue
        tile_path = self.polaris_tp.joinpath(var, self.polaris_stat, depth)
        tile_path.mkdir(parents=True, exist_ok=True)
        urls = [f'{self.polaris_url}/{var}/{self.polaris_stat}/{depth}/{t}.tif' for t in tiles]
        self.downloadTiles(urls, tile_path)


  def processPolaris(self):
    '''
    Merges the POLARIS tiles (or reads the SoilGrids window) of each variable and depth and regrids them onto the
    FastEddy master grid (area average, as for the elevation).
    '''
    tiles = self.polarisTiles() if self.soil_source == 'polaris' else []
    for var in self.polaris_vars:
      for depth in self.polaris_depths:
        mosaic_fp = self.polarisPath(var, depth, 'mosaic')
        if self.soil_source == 'soilgrids':
          if not mosaic_fp.is_file() and not self.polarisPath(var, depth, 'final').is_file():
            self.soilgridsMosaic(var, depth, mosaic_fp)
        elif not mosaic_fp.is_file() and not self.polarisPath(var, depth, 'final').is_file():
          tile_path = self.polaris_tp.joinpath(var, self.polaris_stat, depth)
          sources = [rasterio.open(tile_path.joinpath(f'{t}.tif')) for t in tiles]
          try:
            mosaic, out_trans = merge(sources)
            out_meta = sources[0].meta.copy()
          finally:
            for src in sources:
              src.close()
          out_meta.update({'driver': 'GTiff', 'height': mosaic.shape[1], 'width': mosaic.shape[2],
                           'transform': out_trans, 'compress': 'zstd', 'tiled': True})
          with rasterio.open(mosaic_fp, 'w', **out_meta) as dest:
            dest.write(mosaic)
          logging.info(f'POLARIS mosaic written: {mosaic_fp}')
        self.processRaster(mosaic_fp, self.polarisPath(var, depth, 'final'), Resampling.average)


  def processRaster(self, input_path, output_path, resampling_method):
    '''
    Reprojects, resamples, and clips a raster dataset to exactly match 
    the spatial dimensions and resolution of the FastEddy master grid.

    Parameters
    ----------
    input_path : pathlib.Path
      Path to the source raster.
    output_path : pathlib.Path
      Target path for the fully processed raster.
    resampling_method : rasterio.warp.Resampling
      Method used for pixel interpolation (e.g., bilinear or nearest).
    '''
    xmin, ymin, xmax, ymax = self.db
    
    out_width = int((xmax - xmin) / self.cs)
    out_height = int((ymax - ymin) / self.cs)

    if output_path.is_file():
      with rasterio.open(output_path) as chk:
        if chk.width == out_width and chk.height == out_height:
          logging.info(f'Final data already exists and matches target dimensions at {output_path}.\n')
          return
        else:
          logging.info(f'Stale dimensions found in {output_path}. Overwriting...')
    
    out_transform = transform_from_bounds(xmin, ymin, xmax, ymax, out_width, out_height)
    
    logging.info(f'Processing: {input_path}')
    logging.info(f'Target dimensions: {out_width} x {out_height} pixels')
    
    with rasterio.open(input_path) as src:
      kwargs = src.meta.copy()
      kwargs.update({
        'crs': self.proj_string,
        'transform': out_transform,
        'width': out_width,
        'height': out_height,
        'compress': 'zstd', 
        'tiled': True,
        'bigtiff': 'yes',
        'num_threads': 'all_cpus'
      })

      with rasterio.open(output_path, 'w', **kwargs) as dst:
        for i in range(1, src.count + 1):
          reproject(
            source=rasterio.band(src, i),
            destination=rasterio.band(dst, i),
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=out_transform,
            dst_crs=self.proj_string,
            resampling=resampling_method
          )
          
    logging.info(f'Saved to: {output_path}\n')


  def genCoords(self):
    '''
    Generates the master 1D spatial indices and 2D geographic coordinates 
    (latitude/longitude) arrays for the NetCDF output.
    '''
    logging.info('Creating x and y indices...')
    # Same counts as processRaster (float-safe when the bounds are not integers)
    self.xs = self.db[0] + self.cs * np.arange(int((self.db[2] - self.db[0]) / self.cs))
    self.ys = self.db[1] + self.cs * np.arange(int((self.db[3] - self.db[1]) / self.cs))

    xy = np.meshgrid(self.xs, self.ys)
    logging.info('Index creation complete.')

    logging.info('Transforming to WGS 1984 geodetic coordinates...')
    crs = CRS.from_proj4(self.proj_string)
    t = Transformer.from_crs(crs, 4326) 
    self.lats, self.lons = t.transform(xy[0], xy[1])
    logging.info('Coordinate transform complete.')


  def gisNetCDF(self):
    '''
    Compiles all processed surface fields into a CF-compliant NetCDF-4 dataset,
    ensuring proper array orientations and mapping metadata for FastEddy ingestion.
    '''
    logging.info('Building final NetCDF dataset...')

    with rasterio.open(self.elev_fp) as src:
      elev = src.read(1)
    with rasterio.open(self.nlcd_fp) as src:
      nlcd = src.read(1)
      nlcd_nodata = src.nodata

    elev = np.where(elev == -9999, np.nan, elev)
    elev = fillnodata(elev, mask=~np.isnan(elev))
    # land-cover classes: no-data cells (source no-data value, or 0 outside every CORINE polygon) are filled
    # from neighbouring cells; written as integer class codes
    nlcd_valid = nlcd != nlcd_nodata if nlcd_nodata is not None else np.ones(nlcd.shape, dtype=bool)
    if not np.any(nlcd_valid):
      logging.error(f'The {self.landcover_source} land cover has no data in the domain')
      sys.exit(1)
    nlcd = fillnodata(nlcd, mask=nlcd_valid).astype(np.int16)

    elev = np.flip(elev, axis=0)
    nlcd = np.flip(nlcd, axis=0)

    crs_obj = CRS.from_proj4(self.proj_string)
    wkt_str = crs_obj.to_wkt()
    custom_crs_name = f'{self.dname}_LCC'
    wkt_str = wkt_str.replace('unknown', custom_crs_name, 1)

    data_attrs = {
      'grid_mapping': 'crs',
      'coordinates': 'lat lon'
    }

    data_vars = {
      'x': (['x'], self.xs, {'standard_name': 'projection_x_coordinate', 'units': 'meters'}),
      'y': (['y'], self.ys, {'standard_name': 'projection_y_coordinate', 'units': 'meters'}),
      'cellsize': self.cs,
      
      'lat': (['y', 'x'], self.lats, {'standard_name': 'latitude', 'units': 'degrees_north'}),
      'lon': (['y', 'x'], self.lons, {'standard_name': 'longitude', 'units': 'degrees_east'}),
      
      'elevation': (['y', 'x'], elev, data_attrs),
      'LandCover': (['y', 'x'], nlcd, dict(data_attrs, long_name='land-cover class code',
                                            source=LANDCOVER_SOURCES[self.landcover_source])),
      
      'crs': ([], 0, {
        'grid_mapping_name': 'lambert_conformal_conic',
        'longitude_of_central_meridian': self.proj_lon_0,
        'latitude_of_projection_origin': self.proj_lat_0,
        'standard_parallel': [self.proj_lat_1, self.proj_lat_2],
        'false_easting': 0.0,
        'false_northing': 0.0,
        'crs_wkt': wkt_str,      
        'spatial_ref': wkt_str   
      })
    }

    if self.polaris_url:
      # soil properties per depth interval (soil_depth, y, x); no-data (e.g. open water) filled from neighbours
      units = {'sand': '%', 'silt': '%', 'clay': '%'}
      scale = 0.1 if self.soil_source == 'soilgrids' else 1.0   # SoilGrids texture fractions are in g/kg
      for var in self.polaris_vars:
        layers = []
        for depth in self.polaris_depths:
          with rasterio.open(self.polarisPath(var, depth, 'final')) as src:
            arr = src.read(1).astype(np.float32)
            nodata = src.nodata
          bad = ~np.isfinite(arr) | (arr == -9999)
          if nodata is not None:
            bad |= (arr == nodata)
          arr = np.where(bad, np.nan, arr)
          if np.all(np.isnan(arr)):
            coverage = ' (POLARIS covers CONUS only)' if self.soil_source == 'polaris' else ''
            logging.error(f'{self.soil_source} {var} {depth} cm has no data in the domain{coverage}')
            sys.exit(1)
          arr = fillnodata(arr, mask=~np.isnan(arr)) * scale
          layers.append(np.flip(arr, axis=0))
        data_vars[var] = (['soil_depth', 'y', 'x'], np.stack(layers).astype(np.float32),
                          dict(data_attrs, units=units.get(var, ''),
                               long_name=f'{SOIL_NAMES[self.soil_source]} {var} ({self.polaris_stat}) per depth interval',
                               source=SOIL_SOURCES[self.soil_source]))
      tops = [float(d.split('_')[0]) / 100.0 for d in self.polaris_depths]
      bots = [float(d.split('_')[1]) / 100.0 for d in self.polaris_depths]
      data_vars['soil_depth_top'] = (['soil_depth'], np.array(tops), {'units': 'm', 'long_name': 'top of the depth interval'})
      data_vars['soil_depth_bottom'] = (['soil_depth'], np.array(bots), {'units': 'm', 'long_name': 'bottom of the depth interval'})

    if self.bldgs_fp_url:
      target_raster = self.bldgs_lod1_fp if self.lod == 1 else self.bldgs_fp_fp
      with rasterio.open(target_raster) as src:
        buildings = src.read(1)
      buildings = np.flip(buildings, axis=0)
      data_vars['BuildingHeights'] = (['y', 'x'], buildings, data_attrs)

    ds = xr.Dataset(
      data_vars=data_vars,
      attrs=dict(
        description='FastEddy GIS Dataset',
        Conventions='CF-1.8' 
      )
    )
    ds.to_netcdf(self.nc_fp, format='NETCDF4', engine='netcdf4')

    logging.info('Dataset creation complete. Good bye.')


if __name__ == '__main__':
  if len(sys.argv) == 2:

    param_file = sys.argv[1]
    if Path(param_file).suffix == '.toml':
      with open(param_file, 'rb') as pf:
        setup_params = tomllib.load(pf)
    elif Path(param_file).suffix == '.json':
      with open(param_file, 'r') as pf:
        setup_params = json.load(pf)
      
    domain_name = setup_params['domain_name']
    base_path = Path(setup_params['base_path'])

    base_path.mkdir(parents=True, exist_ok=True)
    
    log_file_name = f'{domain_name}_log.txt'
    log_file = base_path.joinpath(log_file_name)
    
    logging.basicConfig(
      filename=log_file,
      filemode='w',
      level=logging.INFO,
      format='%(asctime)s - %(levelname)s - %(message)s'
    )

    logging.info(f'Preparing FastEddy GIS for {domain_name}')
    FGP = FE_GISPrep(setup_params)
  else:
    print('Usage: python -m FE_GISPrep_Python {path_to_parameter_file}')
    sys.exit(1)

  FGP.bufferedDomain()

  # Check if we already have the elevation data locally before listing the tiles
  if FGP.elev_fp.is_file() or FGP.elev_mp.is_file():
    logging.info('Elevation data already processed or mosaiced. Skipping the tile query.')
  else:
    if FGP.elev_source == 'copernicus':
      elev_urls = FGP.copernicusTiles()
    else:
      elev_urls = FGP.query_tnm_api(FGP.elev_dataset, 'GeoTIFF')
    if elev_urls:
      FGP.downloadTiles(elev_urls, FGP.elev_tp)
      # Store ONLY the specific intersecting files required for this domain
      FGP.elev_files = [FGP.elev_tp.joinpath(url.split('/')[-1]) for url in elev_urls]
    else:
      logging.warning('No elevation tiles found for this domain.')
    

  nlcd_str = str(FGP.nlcd_url).strip()
  if FGP.landcover_source == 'nlcd' and nlcd_str.startswith('http') and nlcd_str.lower().endswith('.zip'):
    FGP.prepareNlcdSource()
    
  if FGP.soil_source == 'polaris':
    FGP.downloadPolaris()

  if FGP.bldgs_fp_url:
    FGP.retrieveBuildingFootprints()
    
  if FGP.lidar_url:
    # Determine which final raster we are trying to build
    if FGP.lod == 1:
      target_bldg_raster = FGP.base_path.joinpath(f'{FGP.dname}_bldgs_lod1_fp.tif')
    else:
      target_bldg_raster = FGP.bldgs_fp_fp

    # If the final raster is already done, skip the API query entirely!
    if target_bldg_raster.is_file():
      logging.info(f'Building data already exists at {target_bldg_raster.name}. Skipping LiDAR API query.')
    else:
      laz_urls = FGP.query_tnm_api('Lidar Point Cloud (LPC)', 'LAZ')
      if laz_urls:
        FGP.downloadTiles(laz_urls, FGP.lidar_tp)
        # Store ONLY the specific intersecting files required for this domain
        FGP.laz_files = [FGP.lidar_tp.joinpath(url.split('/')[-1]) for url in laz_urls]
      else:
        logging.warning('No LiDAR tiles found for this domain.')

  FGP.stitchElev(FGP.elev_mp)
  FGP.bufferedLandCover()
  
  if FGP.lidar_url and FGP.bldgs_fp_url:
    FGP.extractLidarHeights()
  
  FGP.processRaster(FGP.elev_mp, FGP.elev_fp, Resampling.average)
  if FGP.landcover_source == 'corine':
    FGP.rasterizeCorine()
  else:
    FGP.processRaster(FGP.nlcd_bp, FGP.nlcd_fp, FGP.nlcd_resampling)

  if FGP.polaris_url:
    FGP.processPolaris()
  
  if FGP.bldgs_fp_url:
    FGP.processBuildings()

  FGP.genCoords()
  FGP.gisNetCDF()