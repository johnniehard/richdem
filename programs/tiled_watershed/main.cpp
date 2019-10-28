#ifndef _richdem_lindsay2016_hpp_
#define _richdem_lindsay2016_hpp_

#include "richdem/common/Array2D.hpp"
#include "richdem/common/ProgressBar.hpp"
#include "richdem/common/grid_cell.hpp"
#include "richdem/common/logger.hpp"
#include "richdem/common/timer.hpp"
#include "richdem/tiled/A2Array2D.hpp"
#include <limits>
#include <sys/stat.h>
#include <iostream>

using namespace std;
using namespace richdem;

namespace richdem {

template <class elev_t>
void Watershed(A2Array2D<elev_t> &flowdir, A2Array2D<elev_t> &flowacc, int x, int y, int cache_size) {
  cerr << "Starting" << endl;
  
  // Create temporary folder
  mkdir("tmp", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  // Watershed
  mkdir("tmp/watershed", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  A2Array2D<uint32_t> watershed("tmp/watershed/", flowdir.stdTileWidth(), flowdir.stdTileHeight(), flowdir.widthInTiles(), flowdir.heightInTiles(), cache_size);
  watershed.setAll((int)0);

  // Visited
  mkdir("tmp/visited", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  A2Array2D<uint32_t> visited("tmp/visited/", flowdir.stdTileWidth(), flowdir.stdTileHeight(), flowdir.widthInTiles(), flowdir.heightInTiles(), cache_size);
  visited.setAll(false);

  // Start timer
  Timer overall;
  overall.start();

  // queue
  std::queue<std::vector<int>> q;

  std::vector<int> coords(2);
  coords[0] = x;
  coords[1] = y;

  auto flowacc_val = flowacc(x, y);

  cout << "flowacc value at " << x << "," << y << ": " << flowacc_val << endl;

  // add cell to queue
  q.push(coords);

  int cells = 0;
  // while queue is not empty
  while(!q.empty()){
    cells++;

    // get cell from queue
    const auto c = q.front();
    q.pop();

    int x = c[0];
    int y = c[1];

    auto value = flowdir(x, y);

    watershed(x, y) = 1;
    visited(x, y) = true;

    // look at cells neighbours
    for (int n = 1; n <= 8; n++) {
       // continue if cell flows in to n, n will not flow into cell
      if(value == n){
        continue;
      }

      const int nx = x + dx[n];
      const int ny = y + dy[n];

      if(visited(nx, ny)){
        continue;
      }

      auto& n_value = flowdir(nx, ny);

      // check if neighbour flows in to cell
      // if so mark it as part of the watershed and add neighbour to queue
      if(n_value == d8_inverse[n]){
          std::vector<int> coords(2);
          coords[0] = nx;
          coords[1] = ny;
          q.push(coords);
      }
    }
    // cout << "\r" << "cells: " << cells << flush;
  }

  assert(cells == flowacc_val);

  cout << endl;

  // Print final time
  cerr << "Wall-time = " << overall.stop() << endl;
}

} // namespace richdem

int main(int argc, char** argv) {
  if (argc != 5) {
    cerr << "Usage: ./watershed.exe layoutfile_flowdir layoutfile_flowacc x y" << endl;
    cerr << "Layout file should be for flow direction raster" << endl;
    return 1;
  }

  int cache_size = 256;
  A2Array2D<double> flowdir(argv[1], cache_size);
  A2Array2D<double> flowacc(argv[2], cache_size);

  int x = stoi(argv[3]);
  int y = stoi(argv[4]);

  // 736 440

  Watershed(flowdir, flowacc, x, y, cache_size);
}

#endif
