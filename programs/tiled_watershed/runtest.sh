#!/bin/bash

./watershed.exe ../out_flatres/layout.layout ../out_accum/layout.layout $1 $2 | node mangle.js | ogr2ogr -f GeoJSON -s_srs EPSG:3006 -t_srs EPSG:4326 /vsistdout/ /vsistdin/