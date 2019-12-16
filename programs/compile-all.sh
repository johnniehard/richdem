#! /bin/bash

make -C tiler
make -C tiled_lindsay_breaching
make -C tiled_watershed
make -C parallel_priority_flood
make -C tiled_flat_resolution
make -C parallel_d8_accum