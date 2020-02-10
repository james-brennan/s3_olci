"""
s3_pre_process.py
This re-projects a S3_SYN image to MODIS tile and grid
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
import shapely.wkt
import json
from shapely.geometry import shape, GeometryCollection
import glob

if __name__ =="__main__":

    sen3file = sys.argv[1]
    #sen3file = "S3_SYN_20181201T043035_20181201T043237_20181218T164959_0121_038_304_1980"
    sen3filep = Path(sen3file)
    # Figure the MODIS tiles for this image
    with open(sen3filep/"polygon.wkt", 'r') as f:
        poly = f.read()
    olci_swath =  shapely.wkt.loads(poly)
    # get MODIS tile geometries
    mod = '/home/users/jbrennan01/DATA2/BAFMS/data_preprocessing/modis-tiles.geojson'
    with open(mod) as f:
      features = json.load(f)["features"]
    # NOTE: buffer(0) is a trick for fixing scenarios where polygons have overlapping coordinates
    MODIS_TILES = GeometryCollection([shape(feature["geometry"]).buffer(0) for feature in features])
    tiles = np.array([f"h{f['properties']['ih']}v{f['properties']['iv']}" for f in features])

    # ---------- Find intersecting MODIS tiles ----------------------
    intersects = [olci_swath.intersects(tile) for tile in MODIS_TILES]
    idx = np.where(intersects)
    which_tiles = tiles[idx]
    # get amount of intersection

    # --------------- Start processing ------------------------------   
    # create tmp working dir within the file dir
    data_dir = str(sen3filep) + '/'
    work_dir = str(sen3filep) + '/'+ 'tmp'+'/'
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    # 1. Load lon and lat data for preparation of making VRTS
    geo = netCDF4.Dataset(data_dir+"geolocation.nc", 'r')
    lon = geo['lon'][:]
    lat = geo['lat'][:]
    # get dims of raster
    ysize, xsize = lon.shape

    # Save this lon/lat to TIFFs and VRT
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(
        work_dir+"lon.tif",
        xsize,
        ysize,
        1,
        gdal.GDT_Float32, )
    band = dataset.GetRasterBand(1).WriteArray(lon.data)
    dataset = None
    dataset = driver.Create(
        work_dir+"lat.tif",
        xsize,
        ysize,
        1,
        gdal.GDT_Float32, )
    band = dataset.GetRasterBand(1).WriteArray(lat.data)
    dataset = None


    # *-- Construct and save VRTs --*
    lon_vrt = f"""<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">
      <SRS>GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4326"]]</SRS>
  <VRTRasterBand dataType="Float32" band="1">
    <SimpleSource>
      <SourceFilename relativeToVRT="1">lon.tif</SourceFilename>
      <SourceBand>1</SourceBand>
      <SourceProperties RasterXSize="{xsize}" RasterYSize="{ysize}" DataType="Float32" BlockXSize="256" BlockYSize="256" />
      <SrcRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}" />
      <DstRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}" />
    </SimpleSource>
  </VRTRasterBand>
</VRTDataset>"""

    lat_vrt = f"""<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">
      <SRS>GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4326"]]</SRS>
  <VRTRasterBand dataType="Float32" band="1">
    <SimpleSource>
      <SourceFilename relativeToVRT="1">lat.tif</SourceFilename>
      <SourceBand>1</SourceBand>
      <SourceProperties RasterXSize="{xsize}" RasterYSize="{ysize}" DataType="Float32" BlockXSize="256" BlockYSize="256" />
      <SrcRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}" />
      <DstRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}" />
    </SimpleSource>
  </VRTRasterBand>
</VRTDataset>"""

    # save VRT files now
    with open(work_dir+"lon.vrt", "w") as text_file:
        text_file.write(lon_vrt)
    with open(work_dir+"lat.vrt", "w") as text_file:
        text_file.write(lat_vrt)

    # Pull out the data we need and put into a .tif
    dataset = driver.Create(
        work_dir+"data.tif",
        xsize,
        ysize,
        8,
        gdal.GDT_Float32, )

    # *-- Reflectance --*
    #            blue    green   red     NIR
    channels = ["Oa03", "Oa06", "Oa08", "Oa18"]
    for band, chn in enumerate(channels):
        """
        Load data for this band
        """
        ds = netCDF4.Dataset(data_dir+f"Syn_{chn}_reflectance.nc")
        data = ds.variables[f'SDR_{chn}'][:]
        bnd = dataset.GetRasterBand(band+1).WriteArray(data.data)
    """
    *-- QA --*
    Make a QA band from the data
    for now just use the mask provided
    """
    band +=1
    qa = ~data.mask
    bnd = dataset.GetRasterBand(band+1).WriteArray(qa)
    #dataset = None
    """
    *-- Angles --*
    First SZA, VZA, RAA:
    Working on getting angles out of shitty format
    Got to be speed-up for this?
    """
    angs = netCDF4.Dataset(data_dir+"tiepoints_olci.nc", 'r')
    #latS, lonS = angs['OLC_TP_lat'][:], angs['OLC_TP_lon'][:]
    angles = {'OLC_VZA':None, 'SZA':None, 'OLC_VAA':None, 'SAA':None}
    for key in angles.keys():
        an_ = angs[key][:]
        an_ = an_.reshape((-1, ysize)).T
        ang=  resize(an_, (ysize, xsize), order=0, preserve_range=True)
        # check this is the right way around
        angles[key]=ang
    band +=1
    bnd = dataset.GetRasterBand(band+1).WriteArray(angles['OLC_VZA'])
    bnd = dataset.GetRasterBand(band+2).WriteArray(angles['SZA'])
    bnd = dataset.GetRasterBand(band+3).WriteArray(angles['OLC_VAA'] - angles['SAA'])
    # Write to file
    dataset = None
    # Write the main VRT with this data
    band_str = ""
    for band in range(1, 9):
        vrt_raster_band_tmpl = f"""<VRTRasterBand band="{band}" datatype="Float32">
        <SimpleSource>
          <SourceFilename relativeToVRT="1">data.tif</SourceFilename>
          <SourceBand>{band}</SourceBand>
          <SourceProperties RasterXSize="{xsize}" RasterYSize="{ysize}" DataType="Float32" BlockXSize="256" BlockYSize="256" />
          <SrcRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}" />
          <DstRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}" />
        </SimpleSource>
      </VRTRasterBand>
        """
        band_str += vrt_raster_band_tmpl
    # and
    lon_file = work_dir + "lon.tif"
    lat_file = work_dir + "lat.tif"
    vrt_main_tmpl= f"""<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">
       <metadata domain="GEOLOCATION">
         <mdi key="X_DATASET">{lon_file}</mdi>
         <mdi key="X_BAND">1</mdi>
         <mdi key="Y_DATASET">{lat_file}</mdi>
         <mdi key="Y_BAND">1</mdi>
         <mdi key="PIXEL_OFFSET">0</mdi>
         <mdi key="LINE_OFFSET">0</mdi>
         <mdi key="PIXEL_STEP">1</mdi>
         <mdi key="LINE_STEP">1</mdi>
       </metadata>
           {band_str}
</VRTDataset>"""
    # save it
    with open(work_dir+"data.vrt", "w") as text_file:
        text_file.write(vrt_main_tmpl)

    # make a real geotiff now
    g = gdal.Warp(work_dir+"one.tif", work_dir+"data.vrt", dstSRS="EPSG:4326", geoloc=True)
    """
    Perform re-projection to the MODIS grid
    for each tile
    """
    for the_tile in which_tiles:
        # Get the modis reference tile e.g. from MCD64
        t = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_Monthly_500m_DB_BA:Burn Date'
        # form a filename
        keep = sen3filep.stem.split("_")[:8]
        start_str =  "".join([keep[0], "_", keep[3], "_", keep[-1]])
        filename = start_str+f"_{the_tile}.tif"
        outdir = f'/work/scratch/jbrennan01/data/S3/intermediate/{the_tile}/'
        try:
            mcd64dir = f'/home/users/jbrennan01/DATA2/TColBA/input_products/MCD64/{the_tile}/2005/'
            rr = glob.glob(mcd64dir+"*hdf")[0]
            modis_ds = gdal.Open(t % rr)
            geoTransform = modis_ds.GetGeoTransform()
            minx = geoTransform[0]
            maxy = geoTransform[3]
            maxx = minx + geoTransform[1] * modis_ds.RasterXSize
            miny = maxy + geoTransform[5] * modis_ds.RasterYSize
            # Re-project to WGS84 first -- not sure why but seem to need to
            co = 'COMPRESS=DEFLATE', 'INTERLEAVE=BAND', 'PREDICTOR=2'
            g = gdal.Warp(work_dir+f"{the_tile}.tif", work_dir+"one.tif",
                            dstSRS=modis_ds.GetProjection(),  outputBounds=(minx, miny, maxx, maxy),
                                xRes= 463.312719959778804, yRes=-463.312716551443543,
                        creationOptions=co)
            g = None
            """
            Do some cleanup
            1. Delete files
            2. Move output MOD product
            3. rm unzipped file and dir//
            """
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            shutil.move(work_dir+f"{the_tile}.tif", outdir+filename)
            print("success: ", outdir+filename)
        except:
            print("failed:  ", outdir+filename)

    # Now tidy everything up.. eg delete tmp files
    # Gather directory contents
    contents = [os.path.join(work_dir, i) for i in os.listdir(work_dir)]
    # Iterate and remove each item in the appropriate manner
    [os.remove(i) if os.path.isfile(i) or os.path.islink(i) else shutil.rmtree(i) for i in contents]
    os.rmdir(work_dir)



