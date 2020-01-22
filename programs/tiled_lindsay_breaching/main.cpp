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
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;
using namespace richdem;

namespace richdem {

enum LindsayMode {
  COMPLETE_BREACHING,
  SELECTIVE_BREACHING,
  CONSTRAINED_BREACHING
};

enum LindsayCellType { UNVISITED, VISITED, EDGE };

/**
  @brief  Breach and fill depressions (EXPERIMENTAL)
  @author John Lindsay, implementation by Richard Barnes (rbarnes@umn.edu)

    Depression breaching drills a path from a depression's pit cell (its lowest
    point) along the shortest path to the nearest cell outside the depression to
    have the same or lower elevation.

    Several modes are available including:

      *Complete Breaching:    All depressions are entirely breached.
      *Selective Breaching:   Depressions are breached provided the breaching
                              path is not too long nor too deep. That which
                              cannot be breached is filled. Breaching only takes
                              place if the path meets the criteria.
      *Constrained Breaching: A braching path is drilled as long and as deep as
                              permitted, but no more.

    NOTE: It is possible these three modes should be split into different
          functions.

  @param[in,out] &elevations   A grid of cell elevations
  @param[in]     mode          A `LindsayMode` value of `COMPLETE_BREACHING`,
                              `SELECTIVE_BREACHING`, or `CONSTRAINED_BREACHING`.
  @param[in]     eps_gradients If True, then epsilon gradients are applied to
                               breaching paths and depressions to ensure
                               drainage.
  @param[in]     fill_depresssions If True, then depressions are filled.
  @param[in]     maxpathlen    Maximum length of a breaching path
  @param[in]     maxdepth      Maximum depth of a breaching path

  @pre
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM. Note that the _NoData_ value is assumed to
       be a negative number less than any actual data value.

  @return The breached DEM.

  @correctness
    The correctness of this command is determined by inspection and simple unit
    tests.
*/

const uint32_t NO_BACK_LINK = std::numeric_limits<uint32_t>::max();

template <class T>
std::vector<T> read_backwards(std::istream &is, int64_t offset, int size) {
    std::vector<T> buffer;

    buffer.resize(size);

    is.seekg(offset);
    is.read((char*)&buffer[0], size * sizeof(T));
    std::reverse(buffer.begin(), buffer.end());
    return buffer;
}

template <class elev_t>
int elevation2index(elev_t elevation, elev_t minElevation, elev_t maxElevation) {
  int v = (int)((10000 - 1) * (double)(elevation - minElevation)/(double)(maxElevation - minElevation));
  return max(0, min(10000 - 1, v));
}

template <class elev_t>
void Lindsay2016(A2Array2D<elev_t> &dem, int mode, bool eps_gradients,
                 bool fill_depressions, uint32_t maxpathlen, elev_t maxdepth, int cache_size) {
  cerr << "Starting" << endl;
  assert(dem.isReadonly());

  // Move cursor to top left of screen
  std::cout << "\033[0;0;H";
  // Clear screen
  std::cout << "\033[0;J";

  RDLOG_ALG_NAME << "Lindsay2016: Breach/Fill Depressions (EXPERIMENTAL!)";
  RDLOG_CITATION
      << "Lindsay, J.B., 2016. Efficient hybrid breaching-filling sink removal "
         "methods for flow path enforcement in digital elevation models: "
         "Efficient Hybrid Sink Removal Methods for Flow Path Enforcement. "
         "Hydrological Processes 30, 846--857. doi:10.1002/hyp.10648";


  mkdir("tmp", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("tmp/out", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("tmp/elevations", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("tmp/elevations2", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  A2Array2D<float> elevations("tmp/elevations/", dem.stdTileWidth()/4, dem.stdTileHeight()/4, dem.widthInTiles()*4, dem.heightInTiles()*4, cache_size*4*4);
  A2Array2D<float> elevations2("tmp/elevations2/", dem.stdTileWidth()/4, dem.stdTileHeight()/4, dem.widthInTiles()*4, dem.heightInTiles()*4, cache_size*4*4);
  
  for (int y = 0; y < elevations.height(); y++) {
    cout << "\rCopying to native (1 of 2) " << ((int)(100*y/(float)(elevations.height()))) << "%" << flush;
    for (int x = 0;x < elevations.width(); x++) {
      auto elevation = dem(x,y);
      elevations(x, y) = elevation;
    }
  }

  // Clear memory usage
  dem.save_all_tiles();

  cerr << "Determining min and max elevation" << endl;

  auto noDataValue = dem.noData();
  elev_t minElevation = 0;
  elev_t maxElevation = 0;

  for (int y = 0; y < elevations.height(); y++) {
    cout << "\rCopying to native (2 of 2)" << ((int)(100*y/(float)(elevations.height()))) << "%" << flush;
    for (int x = 0;x < elevations.width(); x++) {
      auto elevation = dem(x,y);
      elevations2(x, y) = elevations(x, y) = elevation;

      if (elevation >= 0 && elevation != noDataValue) {
        maxElevation = max(maxElevation, elevation);
        minElevation = min(minElevation, elevation);
      }
    }
  }
  cout << endl;

  cout << "Done" << endl;

  std::vector<uint32_t> flood_array;
  ProgressBar progress;
  Timer overall;

  // Save memory and performance while running breaching
  elevations2.setReadonly(true);

  overall.start();

  uint64_t total_pits = 0;

  // Seed the priority queue

  cout << "Elevation ranges: " << minElevation << "..." << maxElevation << endl; 

  cerr << "Identifying pits and edge cells..." << endl;  

  // List of priority queues. We use multiple ones to reduce the insertion/pop cost (which grows as O(log n))
  vector<GridCellZkB_pq<elev_t>> pqs(elevation2index(maxElevation, minElevation, maxElevation) + 1);
  uint64_t markedAsVisited = 0;

  for (int y = 0; y < dem.height(); y++) {
    cout << "\r" << (int)(100*(double)(y+1)/dem.height()) << "%" << flush;
    //dem.print_cache_debug();

    for (int x = 0; x < dem.width(); x++) {
      auto& elevation = elevations2(x, y);

      if (elevation == noDataValue) // Don't evaluate NoData cells
        continue;

      if (dem.isEdgeCell(x, y)) { // Valid edge cells go on priority-queue
        pqs[elevation2index(elevation, minElevation, maxElevation)].emplace(x, y, elevation, 0);
        // Note: it is important to use elevations instead of dem here
        // They are identical in the beginning, but we will modify elevations while dem
        // stays the same. This is important as the loop below compares elevation values with neighbours.
        elevations(x,y) = noDataValue;
        markedAsVisited++;
        continue;
      }

      // Determine if this is an edge cell, gather information used to determine
      // if it is a pit cell
      elev_t lowest_neighbour = std::numeric_limits<elev_t>::max();
      for (int n = 1; n <= 8; n++) {
        const int nx = x + dx[n];
        const int ny = y + dy[n];

        // No need for an inGrid check here because edge cells are filtered above

        auto nelev = elevations2(nx, ny);
        // Cells which can drain into NoData go on priority-queue as edge cells
        if (nelev == noDataValue) {
          pqs[elevation2index(elevation, minElevation, maxElevation)].emplace(x, y, elevation, 0);
          elevations(x,y) = noDataValue;
          markedAsVisited++;
          goto nextcell; // VELOCIRAPTOR
        }

        // Used for identifying the lowest neighbour
        lowest_neighbour = std::min(nelev, lowest_neighbour);
      }

      // This is a pit cell if it is lower than any of its neighbours. In this
      // case: raise the cell to be just lower than its lowest neighbour. This
      // makes the breaching/tunneling procedures work better.
      // TODO: Include?
      if (elevation < lowest_neighbour) {
        if (eps_gradients)
          elevation = std::nextafter(lowest_neighbour,
                                     std::numeric_limits<elev_t>::lowest());
        else
          elevation = lowest_neighbour;
      }

      // Since depressions might have flat bottoms, we treat flats as pits. Mark
      // flat/pits as such now.
      if (elevation <= lowest_neighbour) {
        //cellInfo(x, y).pit = true;
        total_pits++; // TODO: May not need this
      }

    nextcell:;
    }
  }
  cout << endl;

  // Since dem is read only right now, this will just discard the cache to let us save some memory
  // dem.save_all_tiles();

  // Mark as readonly to save IO bandwidth
  //dem.save_all_tiles();
  assert(dem.isReadonly());

  // The Priority-Flood operation assures that we reach pit cells by passing
  // into depressions over the outlet of minimal elevation on their edge.
  cerr << "Breaching..." << endl;
  uint64_t done = 0;

  std::fstream breach_order;
  breach_order.open("/tmp/breach_order.binary", std::fstream::out | std::fstream::binary);
  
  vector<pair<int64_t, int64_t>> breach_order_buffer;
  for (int pq_index = 0; pq_index < (int)pqs.size();) {
    while (true) {
      auto& pq = pqs[pq_index];
      if (pq.empty()) break;

      done++;
      if ((done % (1024*128)) == 0){
        cout << "\r" << (int)(( done / (double)((int64_t)dem.width()*dem.height())) * 100) << "%. Loops done:" << done << " pq index " << pq_index << "/" << pqs.size() << " pq size: " << pq.size() << " " << markedAsVisited << "/" << ((int64_t)dem.width()*dem.height()) << flush;
      }
      
      const auto c = pq.top();
      int64_t cc = dem.xyToI(c.x, c.y); // Current cell on the path
      pq.pop();
      
      int64_t backlink = c.backlink_dir == 0 ? NO_BACK_LINK : dem.xyToI(c.x + dx[c.backlink_dir], c.y + dy[c.backlink_dir]);
      breach_order_buffer.push_back({cc, backlink});
      if (breach_order_buffer.size() >= 1024*1024) {
        breach_order.write(reinterpret_cast<const char*>(breach_order_buffer.data()),breach_order_buffer.size() * sizeof(pair<int64_t, int64_t>));
        breach_order_buffer.clear();
      }

      // Looks for neighbours which are either unvisited or pits
      for (int n = 1; n <= 8; n++) {
        const int nx = c.x + dx[n];
        const int ny = c.y + dy[n];

        if (!dem.inGrid(nx, ny))
          continue;

        auto& elevation = elevations(nx, ny);

        // This includes checking for visited cells
        if (elevation == noDataValue)
          continue;

        // The neighbour is unvisited. Add it to the queue
        auto backlink_dir = d8_inverse[n];

        int elevIndex = elevation2index(elevation, minElevation, maxElevation);
        pqs[elevIndex].emplace(nx, ny, elevation, backlink_dir);
        if (elevIndex < pq_index) pq_index = elevIndex;

        // Instead of marking the cell as visited, we just set the elevation to an invalid value
        elevation = noDataValue;

        markedAsVisited++;
      }
    }

    assert(pqs[pq_index].empty());
    // Clear the underlaying priority queue storage. Avoids memory leaks.
    pqs[pq_index] = GridCellZkB_pq<elev_t>();
    pq_index++;
  }

  cout << endl;

  // Write remaining items
  breach_order.write(reinterpret_cast<const char*>(breach_order_buffer.data()),breach_order_buffer.size() * sizeof(pair<int64_t, int64_t>));
  breach_order_buffer.clear();
  breach_order.close();

  // Reduce memory usage
  elevations.save_all_tiles();

  cerr << "Applying breaching" << endl;

  auto breach_order_read = ifstream("/tmp/breach_order.binary");
  breach_order_read.seekg(0, std::ios::end);
  uint64_t totalCells = done;
  done = 0;
  while(done < totalCells) {
    uint64_t remaining = totalCells - done;
    uint64_t toRead = min(remaining, 1024 * (uint64_t)1024);
    vector<pair<int64_t, int64_t>> items = read_backwards<pair<int64_t, int64_t>>(breach_order_read, (int64_t)(remaining - toRead) * sizeof(pair<int64_t, int64_t>), (int)toRead);
    done += toRead;

    cout << "\r" << (int)(( done / (float)totalCells) * 100) << "%. Loops done:" << done << flush;

    for (auto item : items) {
      auto cc = item.first;
      auto link = item.second;
      //auto link = cellInfo(cc).backlink;
      if (link != NO_BACK_LINK) {
        auto elev = elevations2(cc);
        // Ensure the backlink's elevation is strictly lower than this cell's elevation
        elevations2(link) = min(elevations2(link), std::nextafter(elev, std::numeric_limits<elev_t>::lowest()));
      }
    }
  }

  breach_order_read.close();

  cerr << "Wall-time = " << overall.stop() << endl;

  // We will write to dem in the loop below. So we mark it as not readonly anymore
  // This is the first time we actually write to the elevation
  // Before this we can treat it as read-only
  dem.cow("tmp/out/");

  for (int y = 0; y < elevations.height(); y++) {
    cout << "\rCopying to different tile size " << ((int)(100*y/(float)(elevations.height()))) << "%" << flush;
    for (int x = 0;x < elevations.width(); x++) {
      dem(x, y) = elevations2(x,y);
    }
  }
}

} // namespace richdem

int main(int argc, char** argv) {
  // Make sure large file support is enabled
  // Otherwise we will may get bugs when we process files larger than about 2GB
  assert(numeric_limits<streamoff>::max() >= numeric_limits<int64_t>::max());

  if (argc != 3) {
    cerr << "Usage: ./breaching layoutfile outpath" << endl;
    cerr << "outpath must include a %f" << endl;
    return 1;
  }

  assert(string(argv[2]).find_last_of("/") != string::npos);
  assert(string(argv[2]).find("%f") != string::npos);

  int cache_size = 16000*4;
  cout << "Reading tile metadata..." << flush;
  A2Array2D<float> dem(argv[1], cache_size);
  cout << " done" << endl;


  Lindsay2016(dem, LindsayMode::COMPLETE_BREACHING, true, true, std::numeric_limits<uint32_t>::max(), std::numeric_limits<float>::max(), cache_size);
  auto foldername = string(argv[2]).substr(0, string(argv[2]).find_last_of("/")).c_str();
  mkdir(foldername, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  cout << "Saving final output" << endl;
  dem.saveGDAL(argv[2]);
}

#endif
