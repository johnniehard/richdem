#!/bin/bash

X=$1
Y=$2
FUNCTION=$3
THRESHOLD=$4

D8_PNTR="/home/johnnie/kod/flodesapp/localdata/geodata/richdem-out/out_flatres/layout.layout"
FLOW_ACC="/home/johnnie/kod/flodesapp/localdata/geodata/richdem-out/out_accum/layout.layout"
NMD="/home/johnnie/kod/flodesapp/localdata/geodata/NMD/tiled_by_richdem/layout.layout"

MANGLE="node mangle.js"
OGR="ogr2ogr -f GeoJSON -s_srs EPSG:3006 -t_srs EPSG:4326 /vsistdout/ /vsistdin/"

# lldb ./watershed.exe -- $D8_PNTR $FLOW_ACC $NMD $X $Y $FUNCTION $THRESHOLD # | $MANGLE | $OGR
./watershed.exe $D8_PNTR $FLOW_ACC $NMD $X $Y $FUNCTION $THRESHOLD # | $MANGLE | $OGR