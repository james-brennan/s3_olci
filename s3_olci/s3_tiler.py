import json

from pathlib import Path

import gdal

import shapely.wkt

from shapely.geometry import shape, GeometryCollection


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
    tile_geom = get_modis_tile_extent("h18v04", "./modis-tiles.geojson" )
    granules, swaths = find_swaths(300, 2019,  "/home/ucfajlg/temp/S3_BRDF/")
    print(granules)
    print(swaths)
    print(get_swaths_overlap_tile(granules, swaths, tile_geom))
    #print(read_s3syn_polyfile(
    #    "/home/ucfajlg/temp/S3_BRDF/2019/300/S3_SYN_20191027T234335_20191027T234527_20191029T102608_0111_051_016_1800"))