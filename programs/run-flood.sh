#! /bin/bash

export GDAL_DISABLE_READDIR_ON_OPEN=1
INPUT_ROADS=$1
set -e

rm -f runflood_log.txt
touch runflood_log.txt

# Tile roads
dt=$(date '+%d/%m/%Y %H:%M:%S')
echo "tiler starting $dt" | tee -a runflood_log.txt
rm -rf out_roads
mkdir out_roads
~/Documents/flodesapp/richdem/programs/tiler/tiler.exe 256 256 $INPUT_ROADS ./out_roads/%n.tif
dt=$(date '+%d/%m/%Y %H:%M:%S')
echo "tiler done $dt" | tee -a runflood_log.txt

# Raise roads
rm -rf out_raised
mkdir out_raised
python ~/Documents/flodesapp/app/dataprep/raise_roads/main.py ./out_breach ./out_tiler ./out_roads ./out_raised
cp ./out_roads/layout.layout ./out_raised

# Run flood
dt=$(date '+%d/%m/%Y %H:%M:%S')
echo "flood starting $dt" | tee -a runflood_log.txt
rm -rf out_flood
mkdir out_flood
mpirun -n 3 ~/Documents/flodesapp/richdem/programs/parallel_priority_flood/parallel_pf.exe many @evict ./out_raised/layout.layout ./out_flood/%f.tif    # 256 256 $INPUT_DEM ./out_tiler/%n.tif
dt=$(date '+%d/%m/%Y %H:%M:%S')
echo "flood done $dt" | tee -a runflood_log.txt

rm -rf out_flood_diff
mkdir out_flood_diff
python ~/Documents/flodesapp/app/dataprep/diff_flood/main.py $PWD/out_breach $PWD/out_flood $PWD/out_flood_diff

echo "all done" | tee -a runflood_log.txt
