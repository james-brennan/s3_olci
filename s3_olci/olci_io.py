
"""
olci_io.py
Assorted code for handling IO of observations...
Given the m month processing idea this will probably be rather
complex. But can be aided by making lots of vrt files for:
"""
import time
import numpy as np
from kernels import *
import scipy.linalg
import glob
import gdal
import datetime
import copy
from dateutil.relativedelta import relativedelta
import os
from rasterio.windows import Window
import rasterio
from affine import *
import scipy.misc
from skimage.transform import rescale, resize



"""
*-- OLCI related IO functions --*
"""
class OLCI_refl(object):
    """
    This is a class to hold and load OLCI reflectance data
    BUT only data produced into the MODIS format by the s3 preprocessor
    """

    def __init__(self, tile, start_date, end_date, xmin, ymin, xmax, ymax):
        self.tile = tile
        self.start_date = start_date
        self.end_date = end_date
        self.xmin = xmin
        self.ymin = ymin
        self.xmax = xmax
        self.ymax = ymax
        self.datadir = f"/work/scratch-nompiio/jbrennan01/data/s3/{tile}/"

    def load_file(self, the_file):
        """
        This loads the file...
        """
        ds = gdal.Open(the_file)
        data = ds.ReadAsArray(
            yoff=int(self.ymin),
            xoff=int(self.xmin),
            ysize=int(self.ymax - self.ymin),
            xsize=int(self.xmax - self.xmin))
        self.ds = ds
        return data

    def loadData(self):
        """
        Loads the viirs refl
        """
        # cba re-writing...
        tile = self.tile
        xmin = self.xmin
        xmax = self.xmax
        ymin = self.ymin
        ymax = self.ymax
        xs = xmax - xmin
        ys = ymax - ymin
        beginning = self.start_date
        ending = self.end_date
        # figure what dates we need
        ndays = (ending - beginning).days
        dates = np.array(
            [beginning + datetime.timedelta(days=x) for x in range(0, ndays)])
        the_files = np.array(glob.glob(self.datadir + "*.tif"))
        the_files.sort()
        file_mtnh = [
            int(f.split("_")[-1].strip(".tif")[4:6]) for f in the_files
        ]
        file_day = [int(f.split("_")[-1].strip(".tif")[6:]) for f in the_files]
        file_years = [
            int(f.split("_")[-1].strip(".tif")[:4]) for f in the_files
        ]
        dates = np.array([
            datetime.datetime(yr, mt, dy)
            for yr, mt, dy in zip(file_years, file_mtnh, file_day)
        ])
        # limit files by dates
        keep = np.logical_and(dates >= self.start_date, dates < self.end_date)
        the_files = the_files[keep]
        dates = dates[keep]
        idx = np.array([(d - self.start_date).days for d in dates])
        """
        make storage for the data
        """
        refl = np.zeros((ndays, 4, ys, xs))
        qa = np.zeros((ndays, ys, xs))
        sza = np.zeros((ndays, ys, xs))
        vza = np.zeros((ndays, ys, xs))
        raa = np.zeros((ndays, ys, xs))
        # loop over the files and load the data
        for f, j in zip(the_files, idx):
            try:
                data = self.load_file(f)
                # we know what each band is quite easily
                _Oa03, _Oa06, _Oa08, _Oa18, _qa, _vza, _sza, _raa = data
                _qa = _qa.astype(np.bool)

                """
                calculate NDVI for rudimentary snow mask
                """
                NDVI = (_Oa18-_Oa08)/(_Oa18+_Oa08)
                snow = NDVI<=0
                _qa = np.logical_and(_qa, ~snow)
                # new-qa as some bands have different qa
                # e.g. B03 seems to have a better qa
                # true is good
                band_b03 = _Oa03<=0
                _qa = np.logical_and(_qa, ~band_b03)
                # sometimes
                # save these...
                refl[j] = (_Oa03, _Oa06, _Oa08, _Oa18)
                refl[j][:, ~_qa]=-1e04
                qa[j] = _qa
                sza[j] = _sza
                vza[j] = _vza
                raa[j] = _raa
            except:
                pass
        """
        compute the kernels for the angles
        """
        kerns = Kernels(
            vza,
            sza, (raa),
            LiType='Sparse',
            doIntegrals=False,
            normalise=True,
            RecipFlag=True,
            RossHS=False,
            MODISSPARSE=True,
            RossType='Thick',
            nbar=0.0)
        kerns.Ross = kerns.Ross.reshape((ndays, ys, xs))
        kerns.Li = kerns.Li.reshape((ndays, ys, xs))
        kerns.Isotropic = kerns.Isotropic.reshape((ndays, ys, xs))
        """
        make everything into a dictionary to be similar across
        sensors
        """
        data = {}
        data['sza'] = sza
        data['vza'] = vza
        data['raa'] = raa
        data['refl'] = refl
        data['qa'] = qa
        data['kernels'] = kerns
        data['date'] = dates
        self.data = data


if __name__ == "__main__":
    if not True:

        tile = 'h20v08'
        start_date = datetime.datetime(2019, 1, 1)
        end_date = datetime.datetime(2019, 3, 31)
        ymin = 0
        ymax = 100
        xmin = 0
        xmax = 100


        o = OLCI_refl(tile, start_date, end_date, xmin, ymin, xmax, ymax)
        o.loadData()

