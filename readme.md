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
- `nlcd_resampling`: rasterio resampling method for land cover (default `'nearest'`; `'mode'` recommended for coarse `cell_size`).
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
550 MB for the default variables and depths. GeoSpec.py (`soil_opt = 1`) turns sand, silt and clay into USDA texture
classes for the LSM.

## Requirements
Python 3.12+