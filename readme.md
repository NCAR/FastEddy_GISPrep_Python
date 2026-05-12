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

## Output
The output GIS NetCDF file matches the specifications required by FastEddy GeoSpec.py.

## Requirements
Python 3.12+