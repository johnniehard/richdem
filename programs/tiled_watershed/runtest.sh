#!/bin/bash

echo $1
echo $2
echo $3
echo $4

# lldb ./watershed.exe -- /home/johnnie/kod/flodesapp/localdata/geodata/richdem-out/out_flatres/layout.layout /home/johnnie/kod/flodesapp/localdata/geodata/richdem-out/out_accum/layout.layout $1 $2
# ./watershed.exe /home/johnnie/kod/flodesapp/localdata/geodata/richdem-out/out_flatres/layout.layout /home/johnnie/kod/flodesapp/localdata/geodata/richdem-out/out_accum/layout.layout $1 $2 # | node mangle.js | ogr2ogr -f GeoJSON -s_srs EPSG:3006 -t_srs EPSG:4326 /vsistdout/ /vsistdin/

./watershed.exe /home/johnnie/kod/flodesapp/localdata/geodata/richdem-out/out_flatres/layout.layout /home/johnnie/kod/flodesapp/localdata/geodata/richdem-out/out_accum/layout.layout $1 $2 $3 $4