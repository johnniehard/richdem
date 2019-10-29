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
#include <vector>

using namespace std;
using namespace richdem;

namespace richdem {

#define PI 3.14159265

std::vector<int> turn(std::vector<int> v, float rad){
    float x = cos(rad) * (float)v[0] - sin(rad) * (float)v[1];
    float y = sin(rad) * (float)v[0] + cos(rad) * (float)v[1];

    return std::vector<int>{(int)x, (int)y};
}

const int n_min = 1;
const int n_max = 8;

int wrapNeighbour(int a){
  return ((a - n_min) % (n_max - n_min) ) + n_min;
}

// Moore neighbourhood boundary tracing
void traceContour(A2Array2D<bool> &watershed, int x, int y){

  // define list of points
  vector<vector<int>> boundary_cells;

  // add starting point to list
  boundary_cells.push_back(vector<int> {x, y});

  // set current cell to x, y
  vector<int> start_cell = {x, y};
  vector<int> cell = start_cell;

  vector<int> start_dir = {0, -1};
  vector<int> dir = start_dir;

  int count = 0;
  while(true){
    if(cell == start_cell && count > 0 && dir == start_dir){
      break;
    }
    count++;

    int cx = cell[0];
    int cy = cell[1];

    if(count > 0){
      dir = turn(dir, -PI);
    }

    for (int i = 1; i <= 8; i++) {

      if(n_diag[i] || i == 1){
        dir = turn(dir, -(PI / 2));
      }

      const int nx = cx - dir[0];
      const int ny = cy - dir[1];


      if(watershed(nx, ny)){
        // cout << "found value" << endl;
        cell[0] = nx;
        cell[1] = ny;
        boundary_cells.push_back(cell);
        break;
      }
    }

    // escape hatch
    if(count > 10000){
      break;
    }

  }

  std::cout << "x, y" << std::endl;
  for(size_t i=0; i < boundary_cells.size(); i++){
    std::cout << boundary_cells[i][0] << ", " << boundary_cells[i][1] << std::endl;
  }

}

template <class elev_t>
void Watershed(A2Array2D<elev_t> &flowdir, A2Array2D<elev_t> &flowacc, int x, int y, int cache_size) {
  // cerr << "Starting" << endl;
  
  // Create temporary folder
  mkdir("tmp", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  // Watershed
  mkdir("tmp/watershed", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  A2Array2D<bool> watershed("tmp/watershed/", flowdir.stdTileWidth(), flowdir.stdTileHeight(), flowdir.widthInTiles(), flowdir.heightInTiles(), cache_size);
  watershed.setAll(false);

  // Visited
  mkdir("tmp/visited", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  A2Array2D<bool> visited("tmp/visited/", flowdir.stdTileWidth(), flowdir.stdTileHeight(), flowdir.widthInTiles(), flowdir.heightInTiles(), cache_size);
  visited.setAll(false);

  // Start timer
  // Timer overall;
  // overall.start();

  // queue
  std::queue<std::vector<int>> q;

  std::vector<int> coords(2);
  coords[0] = x;
  coords[1] = y;

  auto flowacc_val = flowacc(x, y);

  // cout << "flowacc value at " << x << "," << y << ": " << flowacc_val << endl;

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

    watershed(x, y) = true;
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

  // Print final time

  traceContour(watershed, x, y);
  // cerr << "Wall-time = " << overall.stop() << endl;
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
