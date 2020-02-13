# s3_olci

Package for acquiring Sentinel-3 Synergy imagery, mosiacing to a MODIS MOD09-like product and additional BRDF correction code. 

## Searching and downloading Synergy products

To search for and download the synergy products from the Copernicus SciHub use the dload_s3.py script:

```bash
jbrennan01@host492:[~/DATA2/s3_olci/s3_olci]: python dload_s3.py -h
usage: dload_s3.py [-h] [--polygon POLYGON] [--outputdir OUTPUTDIR]
                   [--user USER] [--password PASSWORD]
                   year doy

Download S3 Synergy L2 surface reflectance data.
https://sentinel.esa.int/web/sentinel/user-guides/sentinel-3-synergy/product-
types/level-2-syn

positional arguments:
  year                  Year to download
  doy                   Day of Year to download

optional arguments:
  -h, --help            show this help message and exit
  --polygon POLYGON     A geojson polygon region to limit file downloads to.
                        Otherwise download all SYN files for the day.
  --outputdir OUTPUTDIR
                        Output directory for downloaded files.
  --user USER           Username on scihub.copernicus.eu
  --password PASSWORD   Password for scihub.copernicus.eu
  ```


## Pre-processing to MODIS tiles

Each downloaded S3_SYN product can then be converted to the MODIS sinusoidal grid and tilling system with

```bash
s3_pre_processor.py S3_SYN_file_dir
```
where `S3_SYN_file_dir` is the S3_SYN directory e.g. S3_SYN_20190828T080636_20190828T080936_20190830T040137_0180_048_306_2880


## Producing a MOD09 product




## BRDF correction codes so far
