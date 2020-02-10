import sys
from sentinelsat.sentinel import SentinelAPI, read_geojson, geojson_to_wkt
from datetime import date
import requests
import os
import logging
import datetime
import shapely.wkt
from shapely.geometry.polygon import Polygon
import cartopy.crs as ccrs
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy
"""
This script downloads one day of S3 SYNGERY imagery
"""

auth= ("jbrennan", "burnedarea")
api = SentinelAPI("jbrennan", "burnedarea")

# Some filters
LANDCOVER_KEEP = 5.0
CLOUD_THROW = 95.0
MAX_LAT = 75.0
MIN_LAT = -70.0



import logging
logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
rootLogger = logging.getLogger()
fileHandler = logging.FileHandler("s3_download.log")
fileHandler.setFormatter(logFormatter)
rootLogger.addHandler(fileHandler)
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)


def query(year, doy):
    """
    This does the heavy lifting process
    """
    date = datetime.datetime.strptime(f"{year}{doy}", "%Y%j")
    date0 = date.strftime("%Y%m%d")
    date1 = (date + datetime.timedelta(days=1)).strftime("%Y%m%d")
    products = api.query(area=None, date=(date0, date1),
                         producttype='SY_2_SYN___')
    keep = {}
    """
    Filter products
    For each product we want to filter by landcover
    As the products seem no use
    """
    for p in products.keys():
        try:
            # get exta info eg landcover percentage
            meta = api.get_product_odata(p, full=True)
            # CONDITIONS
            # 1. Landcover greater than 15 -- too high?
            lc = meta['Land Cover Percentage (%)']
            cond1 = lc>LANDCOVER_KEEP
            # 2. Cloud less than 90%?
            cloud = meta['Cloud Cover Percentage (%)']
            cond2 = cloud<CLOUD_THROW
            # 3. Not poles eg below/above
            # this is a bit more challening
            poly = shapely.wkt.loads(meta['footprint'])
            minx, miny, maxx, maxy = poly.bounds
            notPoles = maxy < MAX_LAT and miny> MIN_LAT
            #notPoles = miny> MIN_LAT
            # test all conditions
            passCond = cond1 and cond2 and notPoles
            if passCond:
                keep[p]=meta
        except:
            pass
    return keep


def downloader(product, year, doy):
    """
    Get information about the product
    """
    meta = product
    pid = product['id']
    fname = f"{meta['Filename']}"
    """
    Sort out directory stuff
    """
    base = '/home/users/jbrennan01/DATA2/BAFMS/tmp/SYN/'
    store = base +f'/{year}/{doy}/'
    # mkdir a dir
    if not os.path.exists(store):
        os.makedirs(store)
    fstr= fname.split("____")[1].split("_LN2")[0]
    fdir = store + f'S3_SYN_{fstr}/'
    if not os.path.exists(fdir):
        os.makedirs(fdir)
    # check whether already downloaded
    did_download = check(fdir)
    if not did_download:
        # Download just the files we need
        files_to_get = ['Syn_Oa03_reflectance.nc',
                       'Syn_Oa06_reflectance.nc',
                       'Syn_Oa08_reflectance.nc',
                       'Syn_Oa18_reflectance.nc',
                       'geolocation.nc',
                       'tiepoints_olci.nc']
        for subfilename in files_to_get:
            # format url
            baseurl = (f"https://scihub.copernicus.eu/dhus/odata"+
                        f"/v1/Products(%27{pid}%27)/Nodes(%27{fname}%27)/Nodes(%27{subfilename}%27)/$value")
            #import pdb; pdb.set_trace()
            r = requests.get(baseurl, auth=auth)
            with open(fdir+subfilename, 'wb') as f:
                f.write(r.content)
        # if this has all worked to plan
        # touch ".completed" or something
        with open(fdir+"polygon.wkt", 'w') as f:
            f.write(meta['footprint'])
        with open(fdir+".downloaded", 'w') as f:
            f.write("done")


def make_plot(products, year, doy):
    """
    This plots to show we've downloaded the files we wanted to
    """
    fig = plt.figure(figsize=[5, 3], facecolor="k")
    ax = plt.axes(projection=ccrs.Sinusoidal())
    ax.coastlines(color='white', linewidth=0.5)
    ax.background_patch.set_facecolor('k')
    for p in products.keys():
        # get polygon
        P = shapely.wkt.loads(products[p]['footprint'])
        # check it's been correctly downloaded
        # CODE
        meta = products[p]
        fname = f"{meta['Filename']}"
        base = '/home/users/jbrennan01/DATA2/BAFMS/tmp/SYN/'
        store = base +f'/{year}/{doy}/'
        fstr= fname.split("____")[1].split("_LN2")[0]
        fdir = store + f'S3_SYN_{fstr}/'
        downloaded=check(fdir)
        if downloaded:
            ax.add_geometries([P], crs=ccrs.PlateCarree(), facecolor="none", edgecolor='green', alpha=0.9, linewidth=0.5)
        else:
            ax.add_geometries([P], crs=ccrs.PlateCarree(), facecolor="none", edgecolor='red', alpha=0.9, linewidth=0.5)
    # make the map
    ax.set_title(f"{year} {doy}", color='white')
    plt.savefig(store+"tiles.png", dpi=150, facecolor='k')


def check(fdir):
    """
    check that the file was downloaded
    """
    return os.path.isfile(fdir+'.downloaded')

if __name__ == "__main__":
    # dont forget to change back!
    start = datetime.datetime(2019, 6, 1)
    stop = datetime.datetime(2019, 8, 1)

    while start <= stop:
        year = start.year
        doy = int(start.strftime("%j"))
        logging.info(f"Downloading synergy for {year}{doy}")
        products = query(year, doy)
        for p in list(products.keys()):
            try:
                downloader(products[p], year, doy)
            except:
                print(f"failed for {p}")
        # make a nice plot of the downloaded tiles...
        make_plot(products, year, doy)
        logging.info(f"Downloaded synergy for {year}{doy}")
        start = start + datetime.timedelta(days=1)



