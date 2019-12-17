#! /bin/bash

INPUT_DEM=$1
set -e

rm -f runall_log.txt
rm -f tools_log.txt
touch runall_log.txt
touch tools_log.txt

dt=$(date '+%d/%m/%Y %H:%M:%S')
echo "tiler starting $dt" | tee -a runall_log.txt
rm -rf out_tiler
mkdir out_tiler
~/Documents/flodesapp/richdem/programs/tiler/tiler.exe 256 256 $INPUT_DEM ./out_tiler/%n.tif
dt=$(date '+%d/%m/%Y %H:%M:%S')
echo "tiler done $dt" | tee -a runall_log.txt

dt=$(date '+%d/%m/%Y %H:%M:%S')
echo "breach starting $dt" | tee -a runall_log.txt
rm -rf out_breach
mkdir out_breach
#perf record --call-graph=fp 
~/Documents/flodesapp/richdem/programs/tiled_lindsay_breaching/breaching.exe ./out_tiler/layout.layout ./out_breach/%f.tif
dt=$(date '+%d/%m/%Y %H:%M:%S')
echo "breach done $dt" | tee -a runall_log.txt

dt=$(date '+%d/%m/%Y %H:%M:%S')
echo "flatres starting $dt" | tee -a runall_log.txt
rm -rf out_flatres
mkdir out_flatres
rm -rf out_flatres_tmp
mkdir out_flatres_tmp
~/Documents/flodesapp/richdem/programs/tiled_flat_resolution/parallel_flats.exe ./out_breach/layout.layout 20000 ./out_flatres_tmp/%f.tif ./out_flatres/%f.tif noflip
dt=$(date '+%d/%m/%Y %H:%M:%S')
echo "flatres done $dt" | tee -a runall_log.txt

dt=$(date '+%d/%m/%Y %H:%M:%S')
echo "accum starting $dt" | tee -a runall_log.txt
rm -rf out_accum
mkdir out_accum
# mpirun -n 3 xterm -e lldb ~/Documents/flodesapp/richdem/programs/parallel_d8_accum/parallel_d8_accum.exe many @evict ./out_flatres/layout.layout ./out_accum/%n.tif
mpirun -n 3 ~/Documents/flodesapp/richdem/programs/parallel_d8_accum/parallel_d8_accum.exe many @evict ./out_flatres/layout.layout ./out_accum/%n.tif
dt=$(date '+%d/%m/%Y %H:%M:%S')
echo "accum done $dt" | tee -a runall_log.txt

echo "all done" | tee -a runall_log.txt
