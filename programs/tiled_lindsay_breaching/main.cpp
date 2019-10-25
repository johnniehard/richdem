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
template <class elev_t>
void Lindsay2016(A2Array2D<elev_t> &dem, int mode, bool eps_gradients,
                 bool fill_depressions, uint32_t maxpathlen, elev_t maxdepth) {
  cerr << "Starting" << endl;
  RDLOG_ALG_NAME << "Lindsay2016: Breach/Fill Depressions (EXPERIMENTAL!)";
  RDLOG_CITATION
      << "Lindsay, J.B., 2016. Efficient hybrid breaching-filling sink removal "
         "methods for flow path enforcement in digital elevation models: "
         "Efficient Hybrid Sink Removal Methods for Flow Path Enforcement. "
         "Hydrological Processes 30, 846--857. doi:10.1002/hyp.10648";

  const uint32_t NO_BACK_LINK = std::numeric_limits<uint32_t>::max();

  mkdir("tmp", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("tmp/backlinks", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("tmp/visited", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("tmp/pits", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("tmp/out", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  dem.cow("out/");
  A2Array2D<uint32_t> backlinks("tmp/backlinks/", dem.stdTileWidth(), dem.stdTileHeight(), dem.widthInTiles(), dem.heightInTiles(), 128);
  backlinks.setAll(NO_BACK_LINK);

  A2Array2D<uint8_t> visited("tmp/visited/", dem.stdTileWidth(), dem.stdTileHeight(), dem.widthInTiles(), dem.heightInTiles(), 128);
  visited.setAll(false);

  A2Array2D<uint8_t> pits("tmp/pits/", dem.stdTileWidth(), dem.stdTileHeight(), dem.widthInTiles(), dem.heightInTiles(), 128);
  pits.setAll(false);

  backlinks.copy_metadata_from(dem);
  visited.copy_metadata_from(dem);
  pits.copy_metadata_from(dem);

  std::vector<uint32_t> flood_array;
  GridCellZk_pq<elev_t> pq;
  ProgressBar progress;
  Timer overall;

  overall.start();

  uint32_t total_pits = 0;

  visited.setAll(LindsayCellType::UNVISITED);

  // Seed the priority queue
  cerr << "Identifying pits and edge cells..." << endl;
  
  for (int y = 0; y < dem.height(); y++) {
    cout << "\r " << 100*(float)y/dem.height();
    for (int x = 0; x < dem.width(); x++) {

      if (dem.isNoData(x, y)) // Don't evaluate NoData cells
        continue;

      if (dem.isEdgeCell(x, y)) { // Valid edge cells go on priority-queue
        pq.emplace(x, y, dem(x, y));
        visited(x, y) = LindsayCellType::EDGE;
        continue;
      }

      // Determine if this is an edge cell, gather information used to determine
      // if it is a pit cell
      elev_t lowest_neighbour = std::numeric_limits<elev_t>::max();
      for (int n = 1; n <= 8; n++) {
        const int nx = x + dx[n];
        const int ny = y + dy[n];

        // No need for an inGrid check here because edge cells are filtered
        // above

        // Cells which can drain into NoData go on priority-queue as edge cells
        if (dem.isNoData(nx, ny)) {
          pq.emplace(x, y, dem(x, y));
          visited(x, y) = LindsayCellType::EDGE;
          goto nextcell; // VELOCIRAPTOR
        }

        // Used for identifying the lowest neighbour
        lowest_neighbour = std::min(dem(nx, ny), lowest_neighbour);
      }

      // This is a pit cell if it is lower than any of its neighbours. In this
      // case: raise the cell to be just lower than its lowest neighbour. This
      // makes the breaching/tunneling procedures work better.
      if (dem(x, y) < lowest_neighbour) {
        if (eps_gradients)
          dem(x, y) = std::nextafter(lowest_neighbour,
                                     std::numeric_limits<elev_t>::lowest());
        else
          dem(x, y) = lowest_neighbour;
      }

      // Since depressions might have flat bottoms, we treat flats as pits. Mark
      // flat/pits as such now.
      if (dem(x, y) <= lowest_neighbour) {
        pits(x, y) = true;
        total_pits++; // TODO: May not need this
      }

    nextcell:;
    }
  }
  cout << endl;

  // The Priority-Flood operation assures that we reach pit cells by passing
  // into depressions over the outlet of minimal elevation on their edge.
  cerr << "Breaching..." << endl;
  int total = dem.height() * dem.width();
  int done = 0;
  while (!pq.empty()) {
    done++;
    cout << "\r" << done << "/" << total;

    const auto c = pq.top();
    pq.pop();

    // This cell is a pit: let's consider doing some breaching
    if (pits(c.x, c.y)) {
      // Locate a cell that is lower than the pit cell, or an edge cell
      uint32_t pathlen = 0;
      auto cc = dem.xyToI(c.x, c.y); // Current cell on the path
      elev_t pathdepth =
          std::numeric_limits<elev_t>::lowest(); // Maximum depth found along
                                                 // the path
      elev_t target_height = dem(c.x, c.y); // Depth to which the cell currently
                                            // being considered should be carved

      if (mode == COMPLETE_BREACHING) {
        // Trace path back to a cell low enough for the path to drain into it,
        // or to an edge of the DEM
        while (cc != NO_BACK_LINK && dem(cc) >= target_height) {
          dem(cc) = target_height;
          cc = backlinks(cc); // Follow path back
          if (eps_gradients)
            target_height = std::nextafter(
                target_height,
                std::numeric_limits<
                    elev_t>::lowest()); // Decrease target depth slightly for
                                        // each cell on path to ensure drainage
        }
      } else {
        // Trace path back to a cell low enough for the path to drain into it,
        // or to an edge of the DEM
        while (cc != NO_BACK_LINK && dem(cc) >= target_height) {
          pathdepth = std::max(
              pathdepth,
              (elev_t)(dem(cc) -
                       target_height)); // Figure out deepest breach necessary
                                        // on path //TODO: CHeck this for issues
                                        // with int8_t subtraction overflow
          cc = backlinks(cc); // Follow path back
          if (eps_gradients)
            target_height = std::nextafter(
                target_height,
                std::numeric_limits<
                    elev_t>::lowest()); // Decrease target depth slightly for
                                        // each cell on path to ensure drainage
          pathlen++; // Make path longer
        }

        // Reset current cell address and height to the pit (start of path)
        cc = dem.xyToI(c.x, c.y);
        target_height = dem(c.x, c.y);

        // The path fits within the limits. "Drill, baby, drill."
        if (pathlen <= maxpathlen && pathdepth <= maxdepth) {
          while (cc != NO_BACK_LINK && dem(cc) >= target_height) {
            dem(cc) = target_height;
            cc = backlinks(cc); // Follow path back
            if (eps_gradients)
              target_height = std::nextafter(
                  target_height,
                  std::numeric_limits<
                      elev_t>::lowest()); // Decrease target depth slightly for
                                          // each cell on path to ensure drainage
          }
        } else if (mode == CONSTRAINED_BREACHING) { // TODO: Refine this with
                                                    // regards to the paper
          elev_t current_height = dem(cc);
          while (cc != NO_BACK_LINK && dem(cc) >= target_height) {
            if (pathdepth <= maxdepth)
              dem(cc) = current_height;
            else
              dem(cc) -= pathdepth;
            if (eps_gradients)
              current_height = std::nextafter(
                  current_height, std::numeric_limits<elev_t>::lowest());
            cc = backlinks(cc);
          }
        }
      }

      --total_pits;
      if (total_pits == 0)
        break;
    }

    // Looks for neighbours which are either unvisited or pits
    for (int n = 1; n <= 8; n++) {
      const int nx = c.x + dx[n];
      const int ny = c.y + dy[n];

      if (!dem.inGrid(nx, ny))
        continue;
      if (dem.isNoData(nx, ny))
        continue;
      if (visited(nx, ny) != LindsayCellType::UNVISITED)
        continue;

      const auto my_e = dem(nx, ny);

      // The neighbour is unvisited. Add it to the queue
      pq.emplace(nx, ny, my_e);
      //if (fill_depressions)
      //  flood_array.emplace_back(dem.xyToI(nx, ny));
      visited(nx, ny) = LindsayCellType::VISITED;
      backlinks(nx, ny) = dem.xyToI(c.x, c.y);
    }
  }
  cout << endl;

  // dem.save_all_tiles();
  // visited.save_all_tiles();
  // backlinks.save_all_tiles();
  // pits.save_all_tiles();

  /*
  if (mode != COMPLETE_BREACHING && fill_depressions) {
    RDLOG_PROGRESS << "Flooding...";
    progress.start(dem.numDataCells());
    for (const auto f : flood_array) {
      ++progress;
      auto parent = backlinks(f);
      if (dem(f) <= dem(parent)) {
        if (eps_gradients)
          dem(f) =
              std::nextafter(dem(parent), std::numeric_limits<elev_t>::max());
        else
          dem(f) = dem(parent);
      }
    }
    progress.stop();
  }*/

  cerr << "Wall-time = " << overall.stop() << endl;
}

} // namespace richdem

int main(int argc, char** argv) {
  if (argc != 3) {
    cerr << "Usage: ./a.out layoutfile outpath" << endl;
    cerr << "outpath must include a %f" << endl;
    return 1;
  }

  assert(string(argv[2]).find_last_of("/") != string::npos);
  assert(string(argv[2]).find("%f") != string::npos);

  A2Array2D<double> dem(argv[1], 128);

  Lindsay2016(dem, LindsayMode::COMPLETE_BREACHING, true, true, std::numeric_limits<uint32_t>::max(), std::numeric_limits<double>::max());
  auto foldername = string(argv[2]).substr(0, string(argv[2]).find_last_of("/")).c_str();
  mkdir(foldername, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  dem.saveGDAL(argv[2]);
}

#endif
