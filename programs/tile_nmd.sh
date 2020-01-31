#! /bin/bash

export GDAL_DISABLE_READDIR_ON_OPEN=1
INPUT_DEM=$1
set -e

rm -f tile_nmd.txt
touch tile_nmd.txt

dt=$(date '+%d/%m/%Y %H:%M:%S')
echo "tiler starting $dt" | tee -a tile_nmd.txt
rm -rf out_nmd
mkdir out_nmd
~/Documents/flodesapp/richdem/programs/tiler/tiler.exe 256 256 $INPUT_DEM ./out_nmd/%n.tif
dt=$(date '+%d/%m/%Y %H:%M:%S')
echo "tiler done $dt" | tee -a tile_nmd.txt

echo "all done" | tee -a tile_nmd.txt
