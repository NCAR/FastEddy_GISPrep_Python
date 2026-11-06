# FastEddy GIS Prep - Python

## Overview
The purpose of this code is to prepare a GIS NetCDF file that meets FastEddy input specifications by modifying a simple parameter file and running a single python command.

The workflow is as follows:
1. Download/retrieve all necessary data. If it already exists at the target location, skip this step.
2. Clip to a buffered domain.
3. Reproject, resample, and clip to the target FastEddy LCC projection and domain.
4. Process the resulting .tif files into a single NetCDF dataset.

## Instructions
While in the target directory, clone the remote repository:

git clone https://github.com/NCAR/FastEddy_GISPrep_Python

Navigate to FastEddy_GISPrep_Python. Set up the conda environment by running the following commands:

1. module load conda
2. conda env create -f environment.yml
3. conda activate fasteddy_gisprep

Navigate to src and modify the parameter file to the target domain specifications and set base_path to the target data location. Then run

python -m FE_GISPrep_Python {/path/to/parameter/file}

### Optional parameters
- `proj_lat_1`, `proj_lat_2`, `proj_lat_0`, `proj_lon_0`, `proj_R`: define the Lambert Conformal projection independently of the domain center (`lat_0`, `lon_0`), e.g. to match a WRF parent domain (`TRUELAT1`, `TRUELAT2`, `STAND_LON`, `R = 6370000`). The domain is centered on `lat_0`/`lon_0` with its axes along the projection x/y axes, and the download bounding box is computed from the projected domain.
- `elev_dataset`, `elev_url_filter`: TNM elevation product and URL filter (default: 1 m DEM). For large domains use `'National Elevation Dataset (NED) 1 arc-second'` with `'/Elevation/1/TIFF/'`; only the newest version of each 1x1 degree tile is downloaded.
- `nlcd_resampling` (or `landcover_resampling`): rasterio resampling method for land cover, for every land-cover source (default `'nearest'`; `'mode'` recommended for coarse `cell_size`).
- Set `buildings_fp_url = ''` and `lidar_url = ''` to exclude buildings.

## Output
The output GIS NetCDF file matches the specifications required by FastEddy GeoSpec.py.

## Optional soil data (POLARIS)
For the FastEddy land-surface model, setting `polaris_url` (for example `http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0`)
adds POLARIS soil properties (Chaney et al. 2019, Water Resources Research 55:2916; 30 m, CONUS only). The 1 x 1 degree
tiles covering the buffered domain are downloaded for each variable in `polaris_vars` (default `sand`, `silt`, `clay`, in %),
statistic `polaris_stat` (default `mean`) and depth interval in `polaris_depths` (default all six, `0_5` to `100_200` cm),
merged, and area-averaged onto the FastEddy grid. Each variable is written with dimensions `(soil_depth, y, x)`, together
with `soil_depth_top` and `soil_depth_bottom` (m). Each tile is about 30 MB, so a domain inside one tile needs about
550 MB for the default variables and depths. GeoSpec.py (`lsm_opt = 1`) turns sand, silt and clay into USDA texture
classes for the LSM.

## Data sources outside the US (Europe and global)
The default sources (USGS elevation, NLCD, POLARIS, USGS lidar) cover the US only. Three optional keys select global or
European sources instead. Every one of them is open access and is read directly, with no account. Without these keys the
output is unchanged.

| Key | Values | Source |
|---|---|---|
| `elev_source` | `'tnm'` (default), `'copernicus'` | Copernicus DEM GLO-30: global, 30 m, 1 x 1 degree tiles on AWS open data (`elev_url` overrides the bucket URL). |
| `landcover_source` | `'nlcd'` (default), `'worldcover'`, `'corine'` | ESA WorldCover 2021 v200: global, 10 m, classes 10-100 (`worldcover_url`). CORINE Land Cover 2018: Europe (EEA member and cooperating countries), level-3 codes 111-523 (`corine_url`). |
| `soil_source` | `'polaris'` (implied by `polaris_url`), `'soilgrids'` | ISRIC SoilGrids 2.0: global, 250 m, the same six depth intervals as POLARIS (`soilgrids_url`). |

- **Copernicus DEM:** a surface model, so heights include buildings and tree canopy. Heights are above the EGM2008 geoid.
  Tiles over open ocean do not exist and are skipped.
- **WorldCover:** only the buffered domain is read from the cloud-optimised tiles (HTTP range requests), with no
  full-tile download.
- **CORINE:** the CLC2018 polygons that intersect the buffered domain are queried from the EEA map service and saved to
  `<domain>_corine_buffered.gpkg`. They are then rasterised onto the FastEddy grid: each cell takes the class of the
  polygon containing its centre. The EEA raster download needs an account and is not used.
- **SoilGrids:** sand, silt and clay are read from the ISRIC virtual rasters over the buffered domain, converted from
  g/kg to %, and written exactly as the POLARIS fields. SoilGrids has no data under dense built-up areas and water;
  those cells are filled from their neighbours. `soil_vars`, `soil_stat` and `soil_depths` work like the `polaris_*`
  keys.
- **Land-cover classes:** GeoSpec.py needs a land-cover table with the codes of the chosen source (z0 and, for the
  LSM, the surface properties) and `water_cats` set to its water classes. These are WorldCover 80, and CORINE 511,
  512 and 521-523.
- **Buildings:** Overture footprints are global. USGS lidar is US-only. A LAZ file with no CRS record is assumed to
  be in UTM: NAD83 over North America, ETRS89 over Europe, WGS 84 elsewhere.

Citations and licences:
- Copernicus DEM GLO-30, (c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH 2014-2018, provided under
  COPERNICUS by the European Union and ESA.
- Zanaga et al. (2022), ESA WorldCover 10 m 2021 v200, doi:10.5281/zenodo.7254221, CC-BY 4.0.
- CORINE Land Cover 2018, (c) European Union, Copernicus Land Monitoring Service, European Environment Agency (EEA).
- Poggio et al. (2021), SoilGrids 2.0, SOIL 7:217-240, CC-BY 4.0.

`src/gisprep_murcia.json` is an example for a 20 x 20 km domain over Murcia (Spain) with Copernicus, WorldCover and
SoilGrids.

## Requirements
Python 3.12+