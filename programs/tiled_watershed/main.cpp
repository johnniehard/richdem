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

// Rotate int unit vector by rad radians, i.e. +90deg, -90deg, 180deg
std::vector<int> turn(std::vector<int> v, float rad){
    float x = cos(rad) * (float)v[0] - sin(rad) * (float)v[1];
    float y = sin(rad) * (float)v[0] + cos(rad) * (float)v[1];

    return std::vector<int>{(int)x, (int)y};
}

// Moore neighbourhood boundary tracing
void traceContour(A2Array2D<bool> &watershed, int x, int y, vector<double> transform){

  // define list of points
  vector<vector<int>> boundary_cells;

  // add starting point to list
  boundary_cells.push_back(vector<int> {x, y});

  // set starting point as current cell
  vector<int> start_cell = {x, y};
  vector<int> cell = start_cell;

  vector<int> start_dir;
  // check all neighbours to find empty cell, start towards that cell
  for(int i = 1; i <= 8; i+=2){
    const int nx = x + dx[i];
    const int ny = y + dy[i];
    if(!watershed(nx, ny)){
      start_dir = {dx[i], dy[i]};
    }
  }
  vector<int> dir = start_dir;

  int count = 0;
  while(true){
    // Stopping criteria, stop when reaching start cell for second time
    if(cell == start_cell && count > 0){
      break;
    }
    count++;

    // Backtrack
    if(count > 0){
      dir = turn(dir, PI);
    }

    for (int i = 1; i <= 8; i++) {

      // Turn if on a diagonal or if just backtracked
      if(n_diag[i] || i == 1){
        dir = turn(dir, PI / 2);
      }

      const int nx = cell[0] - dir[0];
      const int ny = cell[1] - dir[1];

      // Check if we've found a value
      // If so set neighbour to active cell and add cell to list of boundary cells
      if(watershed(nx, ny)){
        cell[0] = nx;
        cell[1] = ny;
        boundary_cells.push_back(cell);
        break;
      }
    }

  }
 
  cout << "x, y" << endl;
  for(size_t i=0; i < boundary_cells.size(); i++){

    int raster_x = boundary_cells[i][0];
    int raster_y = boundary_cells[i][1];

    float geo_x = transform[0] + raster_x * transform[1] + raster_y * transform[2];
    float geo_y = transform[3] + raster_x * transform[4] + raster_y * transform[5];

    cout << fixed << setprecision(1) << geo_x << ", " << geo_y << endl;
  }
}

template <class elev_t>
void Watershed(A2Array2D<elev_t> &flowdir, A2Array2D<elev_t> &flowacc, int x, int y, int cache_size) {
 
  // Create temporary folders and datasets
  mkdir("tmp", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  // Watershed
  mkdir("tmp/watershed", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  A2Array2D<bool> watershed("tmp/watershed/", flowdir.stdTileWidth(), flowdir.stdTileHeight(), flowdir.widthInTiles(), flowdir.heightInTiles(), cache_size);
  watershed.setAll(false);
  // Visited
  mkdir("tmp/visited", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  A2Array2D<bool> visited("tmp/visited/", flowdir.stdTileWidth(), flowdir.stdTileHeight(), flowdir.widthInTiles(), flowdir.heightInTiles(), cache_size);
  visited.setAll(false);

  // queue
  std::queue<std::vector<int>> q;

  std::vector<int> coords = {x, y};

  auto flowacc_val = flowacc(x, y);

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
          std::vector<int> coords = {nx, ny};
          q.push(coords);
      }
    }
  }

  assert(cells == flowacc_val);

  traceContour(watershed, x, y, flowacc.getGeotransform());
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

  double geo_x = stod(argv[3]);
  double geo_y = stod(argv[4]);

  // https://gdal.org/api/gdaldataset_cpp.html#_CPPv4N11GDALDataset15GetGeoTransformEPd
  vector<double> transform = flowacc.getGeotransform();
  double raster_x = (1 / (transform[1]*transform[5]-transform[2]*transform[4])) * ( transform[5] * (geo_x - transform[0]) - transform[2] * (geo_y - transform[3]));
  double raster_y = (1 / (transform[1]*transform[5]-transform[2]*transform[4])) * ( transform[1] * (geo_y - transform[3]) - transform[4] * (geo_x - transform[0]));

  Watershed(flowdir, flowacc, raster_x, raster_y, cache_size);
}

#endif
