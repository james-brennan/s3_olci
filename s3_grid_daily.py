"""
s3_grid_daily.py

This produces a daily gridded product from S3A/B images


We can do something a little clever by basically just
listing files and doing this

"""
import gdal
import netCDF4
from skimage.transform import resize
import sys
import zipfile
import os
import numpy as np
import shutil
from pathlib import Path
import glob
import pprint
import datetime 

if __name__ =="__main__":
    _, tile, year, doy = sys.argv
    # makd dir
    outdir = f"/work/scratch-nompiio/jbrennan01/data/s3/{tile}/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # get the choosen date
    choosen_date = datetime.datetime.strptime(f"{year}{doy}", "%Y%j")
    str_date = choosen_date.strftime("%Y%m%d")
    # get files
    ddir = f"/work/scratch/jbrennan01/data/S3/intermediate/{tile}/"
    files = np.array(glob.glob(ddir+"*.tif"))
    dates = np.array([f.split("_")[-3][:8]  for f in files  ])
    idx = np.where(dates==str_date)[0]
    """
    process the files for this_date
    """
    the_files = list(files[idx])
    outfilename = f"S3_SY09_{tile}_{str_date}.tif"
    g = gdal.Warp(outdir+outfilename, the_files, srcNodata=0, creationOptions=["INTERLEAVE=BAND", "TILED=YES", "COMPRESS=DEFLATE", "PREDICTOR=2"])
    g = None
    print(outdir+outfilename, len(the_files))

