from collections import namedtuple
import datetime as dt
import os

from pathlib import Path

from sentinelsat.sentinel import SentinelAPI

import requests

import shapely.wkt
from shapely.geometry.polygon import Polygon
import shutil

DLOAD_OPTS = namedtuple("dload_opts", ["landcover_keep", "cloud_throw",
                        "max_lat", "min_lat"])


def create_outputs(dest_dir, doy, year, fname):
    store = Path(dest_dir) / f"{year}/{doy}"
    store.mkdir(exist_ok=True, parents=True)
    fstr = fname.split("____")[1].split("_LN2")[0]
    fdir = store / f'S3_SYN_{fstr}'
    fdir.mkdir(exist_ok=True, parents=True)
    return fdir


def check_bounds(wkt, max_lat, min_lat):
    minx, miny, maxx, maxy = shapely.wkt.loads(wkt).bounds
    return (maxy < max_lat and miny > min_lat)


def get_s3_fnames(bands):
    """Defines S3 filenames. Takes a list of bands"""
    if type(bands) is not list: bands = [bands]
    files_to_get = ["geolocation.nc", "tiepoints_olci.nc"]
    bandnames = [f"Syn_Oa{band:02}_reflectance.nc" for band in bands]
    files_to_get.extend(bandnames)
    return files_to_get

def download(fname, auth, pid, dest_dir, sel_bands,
             chunks = 1048576):
    for subfilename in get_s3_fnames(sel_bands):
        # format url
        baseurl = (f"https://scihub.copernicus.eu/dhus/odata"+
                    f"/v1/Products(%27{pid}%27)/Nodes(%27{fname}%27)/Nodes(%27{subfilename}%27)/$value")
        r = requests.get(baseurl, auth=auth, stream=True,
                         verify=False)
        if not r.ok:
            raise IOError("Can't start download... [%s]" % source)
        file_size = int(r.headers['content-length'])
        print("Downloading to -> %s" % subfilename)
        print("%d bytes..." % file_size)
        with open(subfilename.replace(".nc", ".partial"), 'wb') as fp:
            cntr = 0
            dload = 0
            for chunk in r.iter_content(chunk_size=chunks):
                if chunk:
                    cntr += 1
                    if cntr > 100:
                        dload += cntr * chunks
                        print("\tWriting %d/%d [%5.2f %%]" %
                                     (dload, file_size,
                                      100. * float(dload) / float(file_size)))
                        cntr = 0

                    fp.write(chunk)
                    fp.flush()
                    os.fsync(fp)
        shutil.move(subfilename.replace("nc", ".partial"), subfilename)


class S3SynergyDowload(object):
    """Sentinel 3 Synergy downloader using SentinelSat
    """
    def __init__(self, username, password, dest_dir,
                 sel_bands = [3,6,8,18], 
                 dload_options=None):

        self.api = SentinelAPI(username, password)
        self.auth = (username, password)
        if dload_options is None:
            self.dload_options = DLOAD_OPTS(5., 95., 75., -70.)
        dest_dir = Path(dest_dir)
        if dest_dir.exists():
            self.dest_dir = dest_dir
        else:
            raise IOError
        self.sel_bands = sel_bands


    def _query(self, doy, year):
        """Uses sentinelsat to query the Query the science hub, or wherever the
        data are kept. Fitlers products by landcover, clouds, etc"""
        date = dt.datetime.strptime(f"{year}{doy}", "%Y%j")
        date0 = date.strftime("%Y%m%d")
        date1 = (date + dt.timedelta(days=1)).strftime("%Y%m%d")
        products = self.api.query(area=None, date=(date0, date1),
                         producttype='SY_2_SYN___')
        selected_products = {k: product for k, product in products.items()
                            if check_bounds(product['footprint'],
                            self.dload_options.max_lat,
                            self.dload_options.min_lat)}
        keep = {}
        for p in selected_products.keys():
            meta = self.api.get_product_odata(p, full=True)
            try:
                # get exta info eg landcover percentage
                # CONDITIONS
                # 1. Landcover greater than 15 -- too high?
                lc = meta["Land Cover Percentage (%)"]
                cond1 = lc > self.dload_options.landcover_keep
                # 2. Cloud less than 90%?
                cloud = meta['Cloud Cover Percentage (%)']
                cond2 = cloud<self.dload_options.cloud_throw
                passCond = cond1 and cond2
                if passCond:
                    keep[p]=meta
            except:
                pass
        print(f"Number of suitable granules {len(keep)}")
        return keep



    def download_data(self, doy, year):
        granules = self._query(doy, year)
        for k, granule in granules.items():
            fname = granule['Filename']
            pid = granule["id"]
            fdir = create_outputs(doy, year, fname)
            download(fname, self.auth, pid, fdir, self.sel_bands)


if __name__ == "__main__":
    s3_dload = S3SynergyDowload("jbrennan", "burnedarea",
            "/home/ucfajlg/temp/S3_BRDF/")
    s3_dload.download_data(300, 2019)