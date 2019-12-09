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

template<class T>
struct MemoryMappedArray2D {
  int width, height;
  T* ptr = nullptr;

  MemoryMappedArray2D() = default;

  MemoryMappedArray2D(string filename, int width, int height) : width(width), height(height) {
    auto numBytes = (size_t)width * (size_t)height * sizeof(T);
    auto fd = open(filename.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    int r = ftruncate(fd, numBytes);
    assert(r == 0);

    ptr = (T*)mmap(nullptr, numBytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    // ptr = (T*)mmap(nullptr, numBytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    // ptr = (T*)malloc(numBytes);
    if (ptr == MAP_FAILED) {
      cout << strerror(errno) << endl;
    }
    assert(ptr != MAP_FAILED);
    close(fd);
    cout << ptr << endl;
  }

  MemoryMappedArray2D(int width, int height) : width(width), height(height) {
    char filename[] = "tmpXXXXXX";
    int fd = mkstemp(filename);
    cout << filename << endl;
    cout << fd << endl;
    unlink(filename);
    auto numBytes = (size_t)width * (size_t)height * sizeof(T);
    int r = ftruncate(fd, numBytes);
    assert(r == 0);

    ptr = (T*)mmap(nullptr, numBytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    // ptr = (T*)mmap(nullptr, numBytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    // ptr = (T*)malloc(numBytes);
    if (ptr == MAP_FAILED) {
      cout << strerror(errno) << endl;
    }
    assert(ptr != MAP_FAILED);
    close(fd);
    cout << ptr << endl;
    //cout << strerror(errno) << endl;
  }

  MemoryMappedArray2D(int width, int height, int fileDescriptor, int fileOffset) : width(width), height(height) {
    auto numBytes = (size_t)width * (size_t)height * sizeof(T);

    ptr = (T*)mmap(nullptr, numBytes, PROT_READ | PROT_WRITE, MAP_SHARED, fileDescriptor, fileOffset);
    // ptr = (T*)mmap(nullptr, numBytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    // ptr = (T*)malloc(numBytes);
    if (ptr == MAP_FAILED) {
      cout << strerror(errno) << endl;
    }
    assert(ptr != MAP_FAILED);
    cout << ptr << endl;
    //cout << strerror(errno) << endl;
  }

  void setAll(T val) {
    for (size_t i = 0; i < (size_t)width * (size_t)height; i++) {
      *(ptr + i) = val;
    }
  }

  ~MemoryMappedArray2D() {
    cout << "Unmapping" << endl;
    munmap(ptr, width * height * sizeof(T));
    //free(ptr);
  }

  T& operator() (int x, int y) {
    assert(x >= 0 && y >= 0 && x < width && y < height);
    return *(ptr + (size_t)y*width + (size_t)x);
  }

  T operator() (int x, int y) const {
    assert(x >= 0 && y >= 0 && x < width && y < height);
    return *(ptr + (size_t)y*width + (size_t)x);
  }

  T& operator() (int64_t i) {
    assert(i >= 0 && i < (int64_t)width*height);
    return *(ptr + (size_t)i);
  }

  T operator() (int64_t i) const {
    assert(i >= 0 && i < (int64_t)width*height);
    return *(ptr + (size_t)i);
  }

  // Delete copy constructor
  // MemoryMappedArray2D(MemoryMappedArray2D& other) = delete;
  MemoryMappedArray2D(const MemoryMappedArray2D&) = delete;
  MemoryMappedArray2D& operator=(const MemoryMappedArray2D&) = delete;
  MemoryMappedArray2D(MemoryMappedArray2D&&) = default;
  MemoryMappedArray2D& operator=(MemoryMappedArray2D&&) = default;
};

template<class T>
struct MemoryMappedGrid {
  int width, height, tileSize, tileSizeLog2;
  int gridWidth, gridHeight;
  vector<MemoryMappedArray2D<T>> grid;

  MemoryMappedGrid(string filename, int width, int height, int tileSize) : width(width), height(height), tileSize(tileSize) {
    assert((tileSize & (tileSize - 1)) == 0 && tileSize > 0);
    gridWidth = (width + tileSize - 1) / tileSize;
    gridHeight = (height+tileSize-1) / tileSize;

    grid.reserve(gridWidth * gridHeight);
    cout << "Creating file" << endl;

    // From example at https://linux.die.net/man/3/open
    auto fd = open(filename.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    size_t numBytes = (size_t)tileSize * (size_t)tileSize * sizeof(T) * (size_t)gridWidth * (size_t)gridHeight;
    assert(ftruncate(fd, numBytes) == 0);

    tileSizeLog2 = 1;
    while((1 << tileSizeLog2) != tileSize) tileSizeLog2++;

    for (int y = 0; y < gridHeight; y++) {
      for (int x = 0; x < gridWidth; x++) {
        cout << x << " " << y << endl;
        grid.emplace_back(tileSize, tileSize, fd, tileSize * tileSize * sizeof(T) * (y * gridWidth + x));
      }
    }

    close(fd);
  }

  void setAll(T val) {
    for (auto& t : grid) t.setAll(val);
  }

  T& operator() (int x, int y) {
    assert(x >= 0 && y >= 0 && x < width && y < height);
    int tx = x >> tileSizeLog2;
    int ty = y >> tileSizeLog2;
    x = x & (tileSize - 1);
    y = y & (tileSize - 1);
    return grid[ty * gridWidth + tx](x, y);
  }

  T operator() (int x, int y) const {
    assert(x >= 0 && y >= 0 && x < width && y < height);
    int tx = x >> tileSizeLog2;
    int ty = y >> tileSizeLog2;
    x = x & (tileSize - 1);
    y = y & (tileSize - 1);
    return grid[ty * gridWidth + tx](x, y);
  }

  T& operator() (int64_t i) {
    return (*this)(i % (int64_t)width, i / (int64_t)width);
  }

  T operator() (int64_t i) const {
    return (*this)(i % (int64_t)width, i / (int64_t)width);
  }
};


const uint32_t NO_BACK_LINK = std::numeric_limits<uint32_t>::max();

#if USE_MMAP
template <typename elev_t>
using BackingArray2D = MemoryMappedArray2D<elev_t>;
#else
template <typename elev_t>
using BackingArray2D = A2Array2D<elev_t>;
#endif

template <class elev_t>
void __attribute__ ((noinline)) traceback(int64_t cc, elev_t target_height, bool eps_gradients, BackingArray2D<elev_t> &elevations, BackingArray2D<uint32_t>& backlinks, BackingArray2D<uint8_t>& visited, BackingArray2D<uint8_t>& pits, BackingArray2D<bool>& endpoint) {
  // Trace path back to a cell low enough for the path to drain into it,
  // or to an edge of the DEM
  // bool isEndpoint;
  while (cc != NO_BACK_LINK && elevations(cc) >= target_height) {
    // if (len == 1) isEndpoint = endpoint((int64_t)cc);
    // if (len == 0) endpoint((int64_t)cc) = true;
    // else endpoint((int64_t)cc) = false;
    //len++;
    elevations(cc) = target_height;
    cc = backlinks(cc); // Follow path back
    if (eps_gradients)
      target_height = std::nextafter(
          target_height,
          std::numeric_limits<
              elev_t>::lowest()); // Decrease target depth slightly for
                                  // each cell on path to ensure drainage
  }

  // if (len > 0) {
  //   len_tot += len && isEndpoint;
  //   len_weight += 1;
  //   hist[min((int)hist.size() - 1, len)] += 1;
  // }
}

template <class T>
std::vector<T> read_backwards(std::istream &is, int size) {
    std::vector<T> buffer;

    buffer.resize(size);

    is.seekg(-size * sizeof(T), std::ios::cur);
    is.read((char*)&buffer[0], size * sizeof(T));
    std::reverse(buffer.begin(), buffer.end());
    return buffer;
}

template <class elev_t>
void Lindsay2016(A2Array2D<elev_t> &dem, int mode, bool eps_gradients,
                 bool fill_depressions, uint32_t maxpathlen, elev_t maxdepth, int cache_size) {
  cerr << "Starting" << endl;


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
  mkdir("tmp/backlinks", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("tmp/visited", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("tmp/pits", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("tmp/out", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("tmp/endpoint", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

#if USE_MMAP
  MemoryMappedArray2D<elev_t> elevations("elevation.binary", dem.width(), dem.height());

  cout << "Do you want to copy the gdal data? [y/n]" << endl;
  string ok;
  cin >> ok;
  if (ok != "n") {
    for (int y = 0; y < elevations.height; y++) {
      cout << "\rCopying to memory mapped file " << ((int)(100*y/(float)(elevations.height))) << "%" << flush;
      for (int x = 0;x < elevations.width; x++) {
        elevations(x, y) = dem(x,y);
      }
    }
    cout << endl;
  }

  cout << "Allocating temporary data" << endl;
  MemoryMappedArray2D<uint32_t> backlinks("backlinks.binary", dem.width(), dem.height());
  MemoryMappedArray2D<uint8_t> visited("visited.binary", dem.width(), dem.height());
  MemoryMappedArray2D<uint8_t> pits("pits.binary", dem.width(), dem.height());
  MemoryMappedArray2D<bool> endpoint("endpoint.binary", dem.width(), dem.height());
#else
  dem.cow("tmp/out/");
  A2Array2D<uint32_t> backlinks("tmp/backlinks/", dem.stdTileWidth(), dem.stdTileHeight(), dem.widthInTiles(), dem.heightInTiles(), cache_size);
  A2Array2D<uint8_t> visited("tmp/visited/", dem.stdTileWidth(), dem.stdTileHeight(), dem.widthInTiles(), dem.heightInTiles(), cache_size);
  A2Array2D<uint8_t> pits("tmp/pits/", dem.stdTileWidth(), dem.stdTileHeight(), dem.widthInTiles(), dem.heightInTiles(), cache_size);
  A2Array2D<bool> endpoint("tmp/endpoint/", dem.stdTileWidth(), dem.stdTileHeight(), dem.widthInTiles(), dem.heightInTiles(), cache_size);
  auto& elevations = dem;
#endif

  endpoint.setAll(false);
  backlinks.setAll(NO_BACK_LINK);
  visited.setAll(LindsayCellType::UNVISITED);
  pits.setAll(false);

  cout << "Done" << endl;

  // backlinks.copy_metadata_from(dem);
  // visited.copy_metadata_from(dem);
  // pits.copy_metadata_from(dem);

  std::vector<uint32_t> flood_array;
  GridCellZk_pq<elev_t> pq;
  ProgressBar progress;
  Timer overall;

  overall.start();

  uint64_t total_pits = 0;

  // Seed the priority queue
  cerr << "Identifying pits and edge cells..." << endl;
  
  vector<int> hist(1000, 0);

  for (int y = 0; y < dem.height(); y++) {
    cout << "\r" << (int)(100*(float)(y+1)/dem.height()) << "%" << flush;
    //dem.print_cache_debug();

    for (int x = 0; x < dem.width(); x++) {
      auto& elevation = elevations(x, y);

      // if (dem.isNoData(x, y)) // Don't evaluate NoData cells
      //   continue;

      if (dem.isEdgeCell(x, y)) { // Valid edge cells go on priority-queue
        pq.emplace(x, y, elevation);
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
        // if (dem.isNoData(nx, ny)) {
        //   pq.emplace(x, y, elevation);
        //   visited(x, y) = LindsayCellType::EDGE;
        //   goto nextcell; // VELOCIRAPTOR
        // }

        // Used for identifying the lowest neighbour
        lowest_neighbour = std::min(elevations(nx, ny), lowest_neighbour);
      }

      // This is a pit cell if it is lower than any of its neighbours. In this
      // case: raise the cell to be just lower than its lowest neighbour. This
      // makes the breaching/tunneling procedures work better.
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
  uint64_t done = 0;
  int64_t pits_left = total_pits;
  double len_tot = 0;
  double len_weight = 0;

  std::fstream breach_order;
  breach_order.open("breach_order.binary", std::fstream::out | std::fstream::binary);
  
  while (!pq.empty()) {
    done++;
    if ((done % (1024*128)) == 0){
      cout << "\r" << (int)(( done / (float)(dem.width()*dem.height())) * 100) << "%. Loops done:" << done << flush;
      //dem.print_cache_debug();
    }
    
    const auto c = pq.top();
    pq.pop();
    breach_order.write(reinterpret_cast<const char*>(&c),sizeof(c)); 

    // This cell is a pit: let's consider doing some breaching
    if (pits(c.x, c.y)) {
      // Locate a cell that is lower than the pit cell, or an edge cell
      //uint32_t pathlen = 0;
      //auto cc = dem.xyToI(c.x, c.y); // Current cell on the path
      //elev_t pathdepth = std::numeric_limits<elev_t>::lowest(); // Maximum depth found along the path
      //elev_t target_height = elevations(c.x, c.y); // Depth to which the cell currently
                                            // being considered should be carved

      if (mode == COMPLETE_BREACHING) {
        //traceback(cc, target_height, eps_gradients, elevations, backlinks, visited, pits, endpoint);
      } else {
        assert(false);
        /*
        // Trace path back to a cell low enough for the path to drain into it,
        // or to an edge of the DEM
        while (cc != NO_BACK_LINK && elevations(cc) >= target_height) {
          pathdepth = std::max(
              pathdepth,
              (elev_t)(elevations(cc) -
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
        target_height = elevations(c.x, c.y);

        // The path fits within the limits. "Drill, baby, drill."
        if (pathlen <= maxpathlen && pathdepth <= maxdepth) {
          while (cc != NO_BACK_LINK && elevations(cc) >= target_height) {
            elevations(cc) = target_height;
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
          elev_t current_height = elevations(cc);
          while (cc != NO_BACK_LINK && elevations(cc) >= target_height) {
            if (pathdepth <= maxdepth)
              elevations(cc) = current_height;
            else
              elevations(cc) -= pathdepth;
            if (eps_gradients)
              current_height = std::nextafter(
                  current_height, std::numeric_limits<elev_t>::lowest());
            cc = backlinks(cc);
          }
       }
       */
      }

      --pits_left;
      if (pits_left == 0){
        cout << "\r" << (int)(( pits_left / (float)total_pits) * 100) << "%. Loops done:" << done << flush;
        break;
      }
    }

    // Looks for neighbours which are either unvisited or pits
    for (int n = 1; n <= 8; n++) {
      const int nx = c.x + dx[n];
      const int ny = c.y + dy[n];

      if (!dem.inGrid(nx, ny))
        continue;
      // if (dem.isNoData(nx, ny))
      //   continue;
      if (visited(nx, ny) != LindsayCellType::UNVISITED)
        continue;

      const auto elevation = elevations(nx, ny);

      // The neighbour is unvisited. Add it to the queue
      pq.emplace(nx, ny, elevation);
      //if (fill_depressions)
      //  flood_array.emplace_back(dem.xyToI(nx, ny));
      visited(nx, ny) = LindsayCellType::VISITED;
      backlinks(nx, ny) = dem.xyToI(c.x, c.y);

      hist[max(0, min((int)hist.size()-1, 500 + (int)(elevation - elevations(c.x, c.y))))] += 1;
    }
  }

  cout << "Lengths: " << (len_tot/len_weight) << " " << len_tot << " " << len_weight << endl;
  cout << endl;
  for (int i = 0; i < hist.size(); i++) {
    cout << i << ": " << hist[i] << endl;
  }

  breach_order.close();

  auto breach_order_read = ifstream("breach_order.binary");
  breach_order_read.seekg(0, std::ios::end);
  uint64_t totalCells = done;
  done = 0;
  while(done < totalCells) {
    uint64_t remaining = totalCells - done;
    uint64_t toRead = min(remaining, 1024 * (uint64_t)1024);
    vector<GridCellZk<elev_t>> items = read_backwards<GridCellZk<elev_t>>(breach_order_read, (int)toRead);
    done += toRead;

    cout << "\r" << (int)(( done / (float)totalCells) * 100) << "%. Loops done:" << done << flush;

    for (auto c : items) {
      auto link = backlinks(c.x, c.y);
      if (link != NO_BACK_LINK) {
        auto elev = elevations(c.x, c.y);
        elevations(link) = min(elevations(link), elev);
      }
    }
  }

  breach_order_read.close();

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
    cerr << "Usage: ./breaching layoutfile outpath" << endl;
    cerr << "outpath must include a %f" << endl;
    return 1;
  }

  MemoryMappedGrid<float> dem_copy("elevation.binary", 1024, 1024, 256);

  // auto arr = MemoryMappedGrid<float>(1024*1024, 100*1024, 4096);
  // arr(500, 500) = 0;
  // for (int y = 0; y < arr.totalHeight; y++) {
  //   for (int x = 0; x < arr.totalWidth; x++) {
  //     arr(x, y) = 1;
  //   }
  // }
  //arr(100* 1024, 1024*1024 - 1) = 10;

  assert(string(argv[2]).find_last_of("/") != string::npos);
  assert(string(argv[2]).find("%f") != string::npos);

  int cache_size = 20000;
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
