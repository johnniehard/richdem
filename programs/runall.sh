#! /bin/bash

INPUT_DEM=$1
set -e

rm -f runall_log.txt
rm -f tools_log.txt
touch runall_log.txt
touch tools_log.txt

dt=$(date '+%Y-%m-%d %H:%M:%S')
echo "tiler starting $dt" >> runall_log.txt
rm -rf out_tiler
mkdir out_tiler
./tiler/tiler.exe 2048 2048 $INPUT_DEM ./out_tiler/%n.tif
dt=$(date '+%Y-%m-%d %H:%M:%S')
echo "tiler done $dt" >> runall_log.txt

dt=$(date '+%Y-%m-%d %H:%M:%S')
echo "breach starting $dt" >> runall_log.txt
rm -rf out_breach
mkdir out_breach
./tiled_lindsay_breaching/breaching.exe ./out_tiler/layout.layout ./out_breach/%f.tif >> tools_log.txt
dt=$(date '+%Y-%m-%d %H:%M:%S')
echo "breach done $dt" >> runall_log.txt

dt=$(date '+%Y-%m-%d %H:%M:%S')
echo "flatres starting $dt" >> runall_log.txt
rm -rf out_flatres
mkdir out_flatres
rm -rf out_flatres_tmp
mkdir out_flatres_tmp
./tiled_flat_resolution/parallel_flats.exe ./out_breach/layout.layout 512 ./out_flatres_tmp/%f.tif ./out_flatres/%f.tif noflip
dt=$(date '+%Y-%m-%d %H:%M:%S')
echo "flatres done $dt" >> runall_log.txt

dt=$(date '+%Y-%m-%d %H:%M:%S')
echo "accum starting $dt" >> runall_log.txt
rm -rf out_accum
mkdir out_accum
mpirun -n 12 ./parallel_d8_accum/parallel_d8_accum.exe many @evict ./out_flatres/layout.layout ./out_accum/%n.tif
dt=$(date '+%Y-%m-%d %H:%M:%S')
echo "accum done $dt" >> runall_log.txt

echo "all done" >> runall_log.txt
