// This program splits a single .tif file into smaller tiles when given a tile size and output directory.
// It also generates a layout file which some other programs in this repository use as input.

// This algorithm is discussed in the manuscript:
//    Barnes, R., 2017. Parallel non-divergent flow accumulation for trillion
//    cell digital elevation models on desktops or clusters. Environmental
//    Modelling & Software 92, 202–212. doi:10.1016/j.envsoft.2017.02.022
#include "gdal_priv.h"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/Layoutfile.hpp"
#include "richdem/common/grid_cell.hpp"
#include "richdem/common/memory.hpp"
#include "richdem/common/timer.hpp"
#include "richdem/common/version.hpp"
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <fstream> //For reading layout files
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream> //Used for parsing the <layout_file>
#include <stack>
#include <string>
#include <vector>

using namespace richdem;
using namespace std;

// We use the cstdint library here to ensure that the program behaves as
// expected across platforms, especially with respect to the expected limits of
// operation for data set sizes and labels. For instance, in C++, a signed
// integer must be at least 16 bits, but not necessarily more. We force a
// minimum of 32 bits as this is, after all, for use with large datasets.
#include <cstdint>

const uint8_t FLIP_VERT = 1;
const uint8_t FLIP_HORZ = 2;

class TileInfo {
private:
  friend class cereal::access;
  template <class Archive> void serialize(Archive &ar) {
    ar(edge, flip, x, y, width, height, gridx, gridy, nullTile, filename,
       outputname, retention, many, analysis);
  }

public:
  uint8_t edge;
  uint8_t flip;
  int32_t x, y, gridx, gridy, width, height;
  bool nullTile;
  bool many;
  std::string filename;
  std::string outputname;
  std::string retention;
  std::string analysis; // Command line command used to invoke everything
  TileInfo() { nullTile = true; }
  TileInfo(std::string filename, std::string outputname, std::string retention,
           int32_t gridx, int32_t gridy, int32_t x, int32_t y, int32_t width,
           int32_t height, bool many, std::string analysis) {
    this->nullTile = false;
    this->edge = 0;
    this->x = x;
    this->y = y;
    this->width = width;
    this->height = height;
    this->gridx = gridx;
    this->gridy = gridy;
    this->filename = filename;
    this->outputname = outputname;
    this->retention = retention;
    this->flip = 0;
    this->many = many;
    this->analysis = analysis;
  }
};

typedef std::vector<std::vector<TileInfo>> TileGrid;

template <class T> void SaveTiles(TileGrid &tiles) {
  const int gridheight = tiles.size();
  const int gridwidth = tiles.front().size();

  cout << endl;
  auto firstTile = tiles[0][0];
  GDALDataset* gdal = (GDALDataset*)GDALOpen(tiles.at(0).at(0).filename.c_str(), GA_ReadOnly);

  for (int y = 0; y < gridheight; y++) {
    for (int x = 0; x < gridwidth; x++) {
      if (tiles[y][x].nullTile)
        continue;

      cout << "\rProcessing " << ((int)(100*(y*gridwidth + x)/(float)(gridwidth*gridheight))) << "% tile " << (y * gridwidth + x) << " of " << (gridwidth * gridheight) << flush;
      auto &tile = tiles.at(y).at(x);

      //auto data = Array2D<T>(tile.filename, false, tile.x, tile.y, tile.width, tile.height);
      auto data = Array2D<T>(gdal, tile.x, tile.y, tile.width, tile.height, false, true);

      // Make sure all tiles are equally large
      auto newData = data;
      newData.resize(firstTile.width, firstTile.height, data.noData());
      for (int y0 = 0; y0 < data.height(); y0++) {
        for (int x0 = 0; x0 < data.width(); x0++) {
          newData(x0, y0) = data(x0, y0);
        }
      }

      newData.saveGDAL(tile.outputname, tile.analysis, tile.x, tile.y);
    }
  }

  cout << endl;
}

// Preparer divides up the input raster file into tiles which can be processed
// independently by the Consumers. Since the tileing may be done on-the-fly or
// rely on preparation the user has done, the Preparer routine knows how to deal
// with both. Once assemebled, the collection of jobs is passed off to Producer,
// which is agnostic as to the original form of the jobs and handles
// communication and solution assembly.
pair<TileGrid, GDALDataType> Preparer(const std::string retention,
                                      const std::string input_file,
                                      const std::string output_name, int bwidth,
                                      int bheight, std::string analysis) {
  Timer timer_overall;
  timer_overall.start();

  TileGrid tiles;
  std::string filename;
  GDALDataType file_type;      // All tiles must have a common file_type
  TileInfo *reptile = nullptr; // Pointer to a representative tile

  std::string output_layout_name = output_name;
  if (output_name.find("%f") != std::string::npos) {
    output_layout_name.replace(output_layout_name.find("%f"), 2, "layout");
  } else if (output_name.find("%n") != std::string::npos) {
    output_layout_name.replace(output_layout_name.find("%n"), 2, "layout");
  } else { // Should never happen
    std::cerr << "E Outputname for mode-many must contain '%f' or '%n'!"
              << std::endl;
    throw std::runtime_error(
        "Outputname for mode-many must contain '%f' or '%n'!");
  }
  LayoutfileWriter lfout(output_layout_name);

  int32_t total_height;
  int32_t total_width;

  // Get the total dimensions of the input file
  try {
    getGDALDimensions(input_file, total_height, total_width, file_type, NULL);
  } catch (...) {
    std::cerr << "E Error getting file information from '" << input_file << "'!"
              << std::endl;
    exit(1);
  }

  // If the user has specified -1, that implies that they want the entire
  // dimension of the raster along the indicated axis to be processed within a
  // single job.
  if (bwidth == -1)
    bwidth = total_width;
  if (bheight == -1)
    bheight = total_height;

  std::cerr << "m Total width =  " << total_width << "\n";
  std::cerr << "m Total height = " << total_height << "\n";
  std::cerr << "m Block width =  " << bwidth << "\n";
  std::cerr << "m Block height = " << bheight << std::endl;
  std::cerr << "m Total cells to be processed = "
            << (total_width * total_height) << std::endl;

  // Create a grid of jobs
  for (int32_t y = 0, gridy = 0; y < total_height; y += bheight, gridy++) {
    tiles.emplace_back(std::vector<TileInfo>());
    lfout.addRow();
    for (int32_t x = 0, gridx = 0; x < total_width; x += bwidth, gridx++) {
      if (retention[0] != '@' && retention.find("%n") == std::string::npos) {
        std::cerr
            << "E In <one> mode '%n' must be present in the retention path."
            << std::endl;
        throw std::invalid_argument("'%n' not found in retention path!");
      }

      if (output_name.find("%n") == std::string::npos) {
        std::cerr << "E In <one> mode '%n' must be present in the output path."
                  << std::endl;
        throw std::invalid_argument("'%n' not found in output path!");
      }

      // Used for '%n' formatting
      std::string coord_string =
          std::to_string(gridx) + "_" + std::to_string(gridy);

      std::string this_retention = retention;
      if (this_retention[0] != '@')
        this_retention.replace(this_retention.find("%n"), 2, coord_string);
      std::string this_output_name = output_name;
      if (this_output_name.find("%n") == std::string::npos) {
        std::cerr << "E Outputname must include '%n' for <one> mode."
                  << std::endl;
        throw std::runtime_error(
            "Outputname must include '%n' for <one> mode.");
      }
      this_output_name.replace(this_output_name.find("%n"), 2, coord_string);

      lfout.addEntry(this_output_name);

      tiles.back().emplace_back(
          input_file, this_output_name, this_retention, gridx, gridy, x, y,
          (total_width - x >= bwidth) ? bwidth : total_width - x, // TODO: Check
          (total_height - y >= bheight) ? bheight : total_height - y, false,
          analysis);
    }
  }

  // If a job is on the edge of the raster, mark it as having this property so
  // that it can be handled with elegance later.
  for (auto &e : tiles.front())
    e.edge |= GRID_TOP;
  for (auto &e : tiles.back())
    e.edge |= GRID_BOTTOM;
  for (size_t y = 0; y < tiles.size(); y++) {
    tiles[y].front().edge |= GRID_LEFT;
    tiles[y].back().edge |= GRID_RIGHT;
  }

  timer_overall.stop();
  std::cerr << "t Preparer time = " << timer_overall.accumulated() << " s"
            << std::endl;

  if (reptile != nullptr) {
    std::cerr << "c Flip horizontal = "
              << ((reptile->flip & FLIP_HORZ) ? "YES" : "NO") << std::endl;
    std::cerr << "c Flip vertical =   "
              << ((reptile->flip & FLIP_VERT) ? "YES" : "NO") << std::endl;
  }
  std::cerr << "c Input data type = " << GDALGetDataTypeName(file_type)
            << std::endl;

  return {tiles, file_type};
}

int main(int argc, char **argv) {
  if (argc != 5) {
    cerr << "Usage ./program <tile width> <tile height> <input file> <output "
            "path>"
         << endl;
    return 1;
  }

  int width = atoi(argv[1]);
  int height = atoi(argv[2]);

  string input = string(argv[3]);

  string output = string(argv[4]);

  assert(width > 0 && height > 0);

  std::string analysis = PrintRichdemHeader(argc, argv);
  TileGrid tiles;
  GDALDataType dataType;
  tie(tiles, dataType) =
      Preparer(output, input, output, width, height, analysis);

  switch (dataType) {
  case GDT_Byte:
    SaveTiles<uint8_t>(tiles);
    break;
  case GDT_UInt16:
    SaveTiles<uint16_t>(tiles);
    break;
  case GDT_Int16:
    SaveTiles<int16_t>(tiles);
    break;
  case GDT_UInt32:
    SaveTiles<uint32_t>(tiles);
    break;
  case GDT_Int32:
    SaveTiles<int32_t>(tiles);
    break;
  case GDT_Float32:
    SaveTiles<float>(tiles);
    break;
  case GDT_Float64:
    SaveTiles<double>(tiles);
    break;
  case GDT_CInt16:
  case GDT_CInt32:
  case GDT_CFloat32:
  case GDT_CFloat64:
    std::cerr << "E Complex types are not supported. Sorry!" << std::endl;
    exit(1);
  case GDT_Unknown:
  default:
    std::cerr << "E Unrecognised data type: " << GDALGetDataTypeName(dataType)
              << std::endl;
    exit(1);
  }
  return 0;
}
