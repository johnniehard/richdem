#! /bin/bash

INPUT_DEM=$1
set -e

echo "tiler"
rm -rf out_tiler
mkdir out_tiler
./tiler/tiler.exe 1024 1024 $INPUT_DEM ./out_tiler/%n.tif

echo "breach"
rm -rf out_breach
mkdir out_breach
./tiled_lindsay_breaching/breaching.exe ./out_tiler/layout.layout ./out_breach/%f.tif

# rm -rf out_flood
# mkdir out_flood
# mpirun -n 4 ./parallel_priority_flood/parallel_pf.exe many @evict ./out_tiler/layout.layout ./out_flood/%n.tif

echo "flatres"
rm -rf out_flatres
mkdir out_flatres
rm -rf out_flatres_tmp
mkdir out_flatres_tmp
./tiled_flat_resolution/parallel_flats.exe ./out_breach/layout.layout 100 ./out_flatres_tmp/%f.tif ./out_flatres/%f.tif noflip

echo "accum"
rm -rf out_accum
mkdir out_accum
mpirun -n 4 ./parallel_d8_accum/parallel_d8_accum.exe many @evict ./out_flatres/layout.layout ./out_accum/%n.tif

# gdal_merge.py out_accum/merged.tif out_accum/*.tif