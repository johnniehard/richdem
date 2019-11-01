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

struct Point {
    int x, y;
    Point(int x, int y) : x(x), y(y) {}
};
bool operator==(const Point& lhs, const Point& rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y;
};

// Rotate int unit vector by rad radians, i.e. +90deg, -90deg, 180deg
Point turn(Point v, float rad){
    float x = cos(rad) * (float)v.x - sin(rad) * (float)v.y;
    float y = sin(rad) * (float)v.x + cos(rad) * (float)v.y;

    return Point((int)x, (int)y);
}

// Moore neighbourhood boundary tracing
void traceContour(A2Array2D<bool> &watershed, Point pour_point, vector<double> transform){

  // define list of points
  vector<Point> boundary_cells;

  // add starting point to list
  boundary_cells.push_back(pour_point);

  // set starting point as current cell
  Point cell = pour_point;

  Point start_dir = Point(0, 0);
  // check all neighbours to find empty cell, start towards that cell
  for(int i = 1; i <= 8; i+=2){
    Point nc = Point(cell.x + dx[i], cell.y + dy[i]);
    if(!watershed(nc.x, nc.y)){
      start_dir = Point(dx[i], dy[i]);
    }
  }
  Point dir = start_dir;

  int count = 0;
  while(true){
    // Stopping criteria, stop when reaching start cell for second time
    if(cell == pour_point && count > 0){
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

      Point nc = Point(cell.x - dir.x, cell.y - dir.y);

      // Check if we've found a value
      // If so set neighbour to active cell and add cell to list of boundary cells
      if(watershed(nc.x, nc.y)){
        cell = nc;
        boundary_cells.push_back(cell);
        break;
      }
    }

  }
 
  cout << "x, y" << endl;
  for(size_t i=0; i < boundary_cells.size(); i++){

    Point bc = boundary_cells[i];

    float geo_x = transform[0] + bc.x * transform[1] + bc.y * transform[2];
    float geo_y = transform[3] + bc.x * transform[4] + bc.y * transform[5];

    cout << fixed << setprecision(1) << geo_x << ", " << geo_y << endl;
  }
}

template <class elev_t>
void Watershed(A2Array2D<elev_t> &flowdir, A2Array2D<elev_t> &flowacc, Point pour_point, int cache_size) {

  assert(pour_point.x > 0 && pour_point.y > 0);
 
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
  std::queue<Point> q;

  // std::vector<int> coords = {x, y};

  auto flowacc_val = flowacc(pour_point.x, pour_point.y);

  // add cell to queue
  q.push(pour_point);

  int cells = 0;
  // while queue is not empty
  while(!q.empty()){
    cells++;

    // get cell from queue
    const auto c = q.front();
    q.pop();

    // int x = c[0];
    // int y = c[1];

    auto value = flowdir(c.x, c.y);

    watershed(c.x, c.y) = true;
    visited(c.x, c.y) = true;

    // look at cells neighbours
    for (int n = 1; n <= 8; n++) {
       // continue if cell flows in to n, n will not flow into cell
      if(value == n){
        continue;
      }

      Point nc = Point(c.x + dx[n], c.y + dy[n]);

      if(visited(nc.x, nc.y)){
        continue;
      }

      auto& n_value = flowdir(nc.x, nc.y);

      // check if neighbour flows in to cell
      // if so mark it as part of the watershed and add neighbour to queue
      if(n_value == d8_inverse[n]){
          q.push(nc);
      }
    }
  }

  assert(cells == flowacc_val);

  traceContour(watershed, pour_point, flowacc.getGeotransform());
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

  Point pour_point = Point((int)raster_x, (int)raster_y);

  Watershed(flowdir, flowacc, pour_point, cache_size);
}

#endif
