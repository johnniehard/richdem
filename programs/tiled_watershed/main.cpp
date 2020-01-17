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
Point operator+(const Point& lhs, const Point& rhs)
{
    return Point(lhs.x + rhs.x, lhs.y + rhs.y);
};
Point operator-(const Point& lhs, const Point& rhs)
{
    return Point(lhs.x - rhs.x, lhs.y - rhs.y);
};

struct Bounds {
  Point min, max;
  Bounds(Point min, Point max) : min(min), max(max) {}
};

// Rotate int unit vector by rad radians, i.e. +90deg, -90deg, 180deg
Point turn(Point v, float rad){
    float x = cos(rad) * (float)v.x - sin(rad) * (float)v.y;
    float y = sin(rad) * (float)v.x + cos(rad) * (float)v.y;

    return Point((int)x, (int)y);
}

const int dx4[4] = { 1, 0, -1, 0 };
const int dy4[4] = { 0, 1, 0, -1 };

// Moore neighbourhood boundary tracing
Bounds traceContour(A2Array2D<bool> &watershed, Point pour_point, vector<double> transform){

  // define list of points
  vector<Point> boundary_cells;

  // set starting point as current cell
  Point cell = pour_point;

  int start_dir = -1;
  // check all neighbours to find empty cell, start towards that cell
  for(int i = 0; i < 4; i++){
    Point nc = Point(cell.x + dx4[i], cell.y + dy4[i]);
    if(!watershed(nc.x, nc.y)){
      start_dir = i;
    }
  }
  assert(start_dir != -1);
  int dir = start_dir;

  // To be able to represent points in-between tiles we double all coordinates
  cell = Point(cell.x*2, cell.y*2);
  cell = cell + Point(dx4[dir], dy4[dir]) + Point(dx4[(dir-1+4)%4], dy4[(dir-1+4)%4]);
  auto start_point = cell;
  dir = (dir+1) % 4;

  int count = 0;
  while(true) {
    if (cell == start_point && count > 0) {
      break;
    }
    count++;

    auto left = (dir + 1) % 4;
    auto right = (dir - 1 + 4) % 4;

    auto front_left = cell + Point(dx4[dir], dy4[dir]) + Point(dx4[left], dy4[left]);
    auto front_right = cell + Point(dx4[dir], dy4[dir]) + Point(dx4[right], dy4[right]);

    auto vleft = watershed(front_left.x/2, front_left.y/2);
    auto vright = watershed(front_right.x/2, front_right.y/2);

    if (vleft && !vright) {
      // Move forwards
      boundary_cells.push_back(cell);
      cell = cell + Point(dx4[dir]*2, dy4[dir]*2);
    } else if (vright) {
      dir = (dir - 1 + 4) % 4;
    } else {
      dir = (dir + 1) % 4;
    }
  }
  /*
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

  }*/

  Bounds bounds = Bounds(Point(INT_MAX, INT_MAX), Point(INT_MIN, INT_MIN));
  
  // Compensate for the doubling of coordinates mentioned above
  transform[1] *= 0.5f;
  transform[4] *= 0.5f;
  transform[2] *= 0.5f;
  transform[5] *= 0.5f;
  
  cout << "x, y" << endl;
  for(size_t i=0; i < boundary_cells.size(); i++){

    Point bc = boundary_cells[i];

    if(bc.x < bounds.min.x){
      bounds.min.x = bc.x;
    }
    if(bc.y < bounds.min.y){
      bounds.min.y = bc.y;
    }
    if(bc.x > bounds.max.x){
      bounds.max.x = bc.x;
    }
    if(bc.y > bounds.max.y){
      bounds.max.y = bc.y;
    }

    float geo_x = transform[0] + bc.x * transform[1] + bc.y * transform[2];
    float geo_y = transform[3] + bc.x * transform[4] + bc.y * transform[5];

    cout << fixed << setprecision(1) << geo_x << ", " << geo_y << endl;
  }

  return bounds;
}

template <class elev_t>
A2Array2D<bool> Watershed(A2Array2D<elev_t> &flowdir, A2Array2D<elev_t> &flowacc, Point pour_point, int cache_size) {

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

  return watershed;
}

template <class elev_t>
void SnapToFlowacc(A2Array2D<elev_t> &flowdir, A2Array2D<elev_t> &flowacc, Point pour_point, int cache_size, int threshold) {

  Point c = pour_point;

  // Follow flow direction until meeting a flow accumulation >= threshold
  while(true){

    int flowacc_val = flowacc(c.x, c.y);

    // Check current cell value against threshold
    if(flowacc_val >= threshold){
      break;
    }

    // look up flow direction
    int dir = flowdir(c.x, c.y);

    // set current cell to downstream neighbour
    c = Point(c.x + dx[dir], c.y + dy[dir]);
  }

  // Print out geo coordinates of cell
  cout << "x, y" << endl;

  vector<double> transform = flowacc.getGeotransform();

  float geo_x = transform[0] + c.x * transform[1] + c.y * transform[2];
  float geo_y = transform[3] + c.x * transform[4] + c.y * transform[5];

  cout << fixed << setprecision(1) << geo_x << ", " << geo_y << endl;
}

void zonalStats(A2Array2D<bool> &watershed, A2Array2D<int> &nmd, Bounds watershed_bounds, vector<double> ws_transform) {

  vector<double> nmd_transform = nmd.getGeotransform();

  int stats[256] = { 0 };

  for (int y = watershed_bounds.min.y; y < watershed_bounds.max.y; y++) {
    // cout << "\r" << (int)(100*(float)(y+1)/watershed_bounds.max.y) << "%" << flush;
    for (int x = watershed_bounds.min.x; x < watershed_bounds.max.x; x++) {
      bool val = watershed(x, y);
      if(val){
        // Get SWEREF coord for watershed xy
        float geo_x = ws_transform[0] + x * ws_transform[1] + y * ws_transform[2];
        float geo_y = ws_transform[3] + x * ws_transform[4] + y * ws_transform[5];
        
        // Get NMD xy for SWEREF coord
        double raster_x = (1 / (nmd_transform[1]*nmd_transform[5]-nmd_transform[2]*nmd_transform[4])) * ( nmd_transform[5] * (geo_x - nmd_transform[0]) - nmd_transform[2] * (geo_y - nmd_transform[3]));
        double raster_y = (1 / (nmd_transform[1]*nmd_transform[5]-nmd_transform[2]*nmd_transform[4])) * ( nmd_transform[1] * (geo_y - nmd_transform[3]) - nmd_transform[4] * (geo_x - nmd_transform[0]));

        // Get NMD value
        int nmd_value = nmd(raster_x, raster_y);

        // Update stats array
        stats[nmd_value] += 1;
      }
    }
  }

  cout << "nmd_klass, count" << endl;
  for (int i = 0; i < 256; i++){
    if(stats[i] > 0){
      cout << i << ", " << stats[i] << endl;
    }
  }
}

} // namespace richdem

void usage(){
  cerr << "Usage: ./watershed.exe layoutfile_flowdir layoutfile_flowacc layoutfile_NMD x y function <threshold>" << endl;
  cerr << "Layout file should be for flow direction raster" << endl;
  cerr << "function can be either 'watershed' or 'snap'" << endl;
  cerr << "pass a threshold value if using snap" << endl;
}

int main(int argc, char** argv) {

  string function = string(argv[6]);

  if(function == "watershed"){
    if (argc != 7) {
      usage();
      return 1;
    }
  }
  else if(function == "snap"){
    if (argc != 8) {
      usage();
      return 1;
    }
  } else {
    usage();
    return 1;
  }

  int cache_size = 256;
  A2Array2D<double> flowdir(argv[1], cache_size);
  A2Array2D<double> flowacc(argv[2], cache_size);

  double geo_x = stod(argv[4]);
  double geo_y = stod(argv[5]);


  // https://gdal.org/api/gdaldataset_cpp.html#_CPPv4N11GDALDataset15GetGeoTransformEPd
  vector<double> transform = flowacc.getGeotransform();
  double raster_x = (1 / (transform[1]*transform[5]-transform[2]*transform[4])) * ( transform[5] * (geo_x - transform[0]) - transform[2] * (geo_y - transform[3]));
  double raster_y = (1 / (transform[1]*transform[5]-transform[2]*transform[4])) * ( transform[1] * (geo_y - transform[3]) - transform[4] * (geo_x - transform[0]));

  Point pour_point = Point((int)raster_x, (int)raster_y);

  if(function == "watershed"){
    A2Array2D<bool> watershed = Watershed(flowdir, flowacc, pour_point, cache_size);
    Bounds watershed_bounds = traceContour(watershed, pour_point, flowacc.getGeotransform());
    // cout << watershed_bounds.min.x <<  " " << watershed_bounds.min.y <<  " " << watershed_bounds.max.x <<  " " << watershed_bounds.max.y << endl;
    cout << "---" << endl;
    A2Array2D<int> nmd(argv[3], cache_size);
    zonalStats(watershed, nmd, watershed_bounds, transform);
  }
  else if(function == "snap"){
    int theshold = stoi(argv[7]);
    SnapToFlowacc(flowdir, flowacc, pour_point, cache_size, theshold);
  }
}

#endif
