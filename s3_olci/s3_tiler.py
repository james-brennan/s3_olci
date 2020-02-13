import json

from pathlib import Path

import gdal

import shapely.wkt

from shapely.geometry import shape, GeometryCollection


def create_vrts(fname, xsize, ysize):
    vrt = f"""<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">
          <SRS>GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4326"]]</SRS>
      <VRTRasterBand dataType="Float32" band="1">
        <SimpleSource>
          <SourceFilename relativeToVRT="1">{fname}</SourceFilename>
          <SourceBand>1</SourceBand>
          <SourceProperties RasterXSize="{xsize}" RasterYSize="{ysize}" DataType="Float32" BlockXSize="256" BlockYSize="256" />
          <SrcRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}" />
          <DstRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}" />
        </SimpleSource>
      </VRTRasterBand>
    </VRTDataset>"""

    def convert_geoloc_tiff(data_dir):
        nc_file = data_dir / "geolocation.nc"
        for dimension in ["lat", "lon"]:
            fname = data_dir/f"{dimension}.tif"
            g = gdal.Warp( fname.as_posix(),
                        f"'NETCDF:{nc_file}:{dimension}'",
                        format="GeoTIFF",
                        output_type=gdal.GDT_Float32)
            xsize = g.RasterXSize
            ysize = g.RasterYSize
            vrt = create_vrts(f"{dimension}.tif", xsize, ysize)
            fname.write_text(vrt)


def read_s3syn_polyfile(s3_granule):
    """Reads the S3 WKT polygon. useful for figuring stuff out"""
    s3_granule = Path(s3_granule)
    return shapely.wkt.loads(s3_granule.read_text())

def get_modis_tile_extent(selected_tile,  modis_tile_list):
    modis_tile_list = Path(modis_tile_list)
    features = json.loads(modis_tile_list.read_text())['features']
    compare = lambda f: f"h{f['properties']['ih']}v{f['properties']['iv']}"
    modis_tile = GeometryCollection([shape(feat['geometry']).buffer(0)
                                    for feat in features
                                    if compare(feat) == selected_tile])
    return modis_tile

def find_swaths(doy, year, data_dir):
    path = Path(data_dir) / f"{year}/{doy}"
    swaths = [read_s3syn_polyfile(f) for f in path.rglob("**/polygon.wkt")]
    granules = [f.parent for f in path.rglob("**/polygon.wkt")]
    return granules, swaths

def get_swaths_overlap_tile(swaths, swath_geoms, tile_geom):
    sel_swaths = [swaths[i] for i,swath in enumerate(swath_geoms)
                  if tile_geom.intersects(swath)]
    return sel_swaths


if __name__ == "__main__":
    h = 18
    v = 2
    # h18v02 seems to be present
    tile_geom = get_modis_tile_extent(f"h{h:02}v{v:02}", "./modis-tiles.geojson" )
    granules, swaths = find_swaths(300, 2019,  "/home/ucfajlg/temp/S3_BRDF/")
    overlaps = get_swaths_overlap_tile(granules, swaths, tile_geom)
    if overlaps:
        print(f"h{h:02}v{v:02} => ", overlaps)
    #print(read_s3syn_polyfile(
    #    "/home/ucfajlg/temp/S3_BRDF/2019/300/S3_SYN_20191027T234335_20191027T234527_20191029T102608_0111_051_016_1800"))