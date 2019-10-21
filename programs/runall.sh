#! /bin/bash

INPUT_DEM=$1

rm -rf out_tiler
mkdir out_tiler
./tiler/tiler.exe 2048 2048 $INPUT_DEM ./out_tiler/%n.tif

rm -rf out_flood
mkdir out_flood
mpirun -n 4 ./parallel_priority_flood/parallel_pf.exe many @evict ./out_tiler/layout.layout ./out_flood/%n.tif

rm -rf out_flatres
mkdir out_flatres
rm -rf out_flatres_tmp
mkdir out_flatres_tmp
./tiled_flat_resolution/parallel_flats.exe ./out_flood/layout.layout 100 ./out_flatres_tmp/%f.tif ./out_flatres/%f.tif noflip

rm -rf out_accum
mkdir out_accum
mpirun -n 4 ./parallel_d8_accum/parallel_d8_accum.exe many @evict ./out_flatres/layout.layout ./out_accum/%n.tif

# gdal_merge.py out_accum/merged.tif out_accum/*.tif