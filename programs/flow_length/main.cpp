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
#include "external_sort/external_sort.hpp"

using namespace std;
using namespace richdem;


//Valid flowdirs are in the range 0-8, inclusive. 0 is the center and 1-8,
//inclusive, are the 8 cells around the central cell. An extra value is needed
//to indicate NoData. Therefore, uint8_t is appropriate.
typedef uint8_t flowdir_t;

struct Node {
    int x, y;
    float height;
    flowdir_t dir;

    bool operator< (const Node& other) const {
        return other.height < height;
    }
};

A2Array2D<int> calculate_flow(A2Array2D<flowdir_t> &directions, int cache_size) {
    /*std::fstream node_order;
    node_order.open("node_order.binary", std::fstream::out | std::fstream::binary);
    vector<Node> buffer;

    for (int y = 0; y < dem.height(); y++) {
        cout << "\r" << (int)(100*(float)(y+1)/dem.height()) << "%" << flush;

        for (int x = 0; x < dem.width(); x++) {
            auto& elevation = dem(x, y);
            auto& dir = directions(x, y);
            if (dir == 0) continue;
            if (!isfinite(elevation)) continue;

            Node node;
            node.x = x;
            node.y = y;
            node.height = elevation;//y*dem.width() + x;
            node.dir = dir;
            buffer.push_back(node);
            
            if (buffer.size() >= 1024*1024) {
                node_order.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(Node));
                buffer.clear();
            }
        }
    }

    node_order.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(Node));
    buffer.clear();
    node_order.close();

    // set split and merge parameters
    external_sort::SplitParams sp;
    external_sort::MergeParams mp;
    sp.mem.size = 10;
    sp.mem.unit = external_sort::MB;
    mp.mem = sp.mem;
    sp.spl.ifile = "node_order.binary";
    mp.mrg.ofile = "node_order_sorted.binary";

    using ValueType = Node;

    // run external sort
    external_sort::sort<ValueType>(sp, mp);

    if (sp.err.none && mp.err.none) {
        std::cout << "File sorted successfully!" << std::endl;
    } else {
        std::cout << "External sort failed!" << std::endl;
        if (sp.err) {
            std::cout << "Split failed: " << sp.err.msg() << std::endl;
        } else {
            std::cout << "Merge failed: " << mp.err.msg() << std::endl;
        }
    }*/

    // node_order.open("node_order.binary", std::fstream::out | std::fstream::binary);
    // vector<Node> buffer;

    A2Array2D<int> in_count("tmp/in_count/", directions.stdTileWidth(), directions.stdTileHeight(), directions.widthInTiles(), directions.heightInTiles(), cache_size);
    in_count.setAll(0);

    for (int y = 0; y < directions.height(); y++) {
        cout << "\r" << (int)(100*(float)(y+1)/directions.height()) << "%" << flush;

        for (int x = 0; x < directions.width(); x++) {
            auto& dir = directions(x, y);
            if (dir == 0) continue;
            if (in_count.in_grid(x + dx[dir], y + dy[dir])) {
                in_count(x + dx[dir], y + dy[dir])++;
            }
        }
    }

    queue<tuple<int,int>> que;
    for (int y = 0; y < directions.height(); y++) {
        cout << "\r" << (int)(100*(float)(y+1)/directions.height()) << "%" << flush;

        for (int x = 0; x < directions.width(); x++) {
            if (in_count(x, y) == 0) que.push({x, y});
        }
    }

    A2Array2D<int> flow_lengths("tmp/flow_length/", directions.stdTileWidth(), directions.stdTileHeight(), directions.widthInTiles(), directions.heightInTiles(), cache_size);
    flow_lengths.copy_metadata_from(directions);
    flow_lengths.setAll(1);

    int record = 0;
    while(!que.empty()) {
        int x, y;
        tie(x, y) = que.front();
        que.pop();
        auto dir = directions(x, y);
        auto nx = x + dx[dir];
        auto ny = y + dy[dir];
        if (!flow_lengths.in_grid(nx, ny)) continue;

        auto current = flow_lengths(x, y);
        auto& next = flow_lengths(nx, ny);
        next = max(next, current + 1);
        record = max(next, record);
        auto& next_in = in_count(nx, ny);
        next_in--;
        if (next_in == 0) que.push({nx, ny});
    }

    /*
    auto node_order_sorted = ifstream("node_order_sorted.binary", std::fstream::binary);
    std::vector<Node> node_read_buffer(1024*1024);

    dem.save_all_tiles();
    directions.save_all_tiles();

    A2Array2D<int> flow_lengths("tmp/flow_length/", dem.stdTileWidth(), dem.stdTileHeight(), dem.widthInTiles(), dem.heightInTiles(), cache_size);
    flow_lengths.copy_metadata_from(dem);
    flow_lengths.setAll(1);

    cout << " Total size " << dem.width() << " " << dem.height() << endl;
    int record = 0;
    do {
        node_order_sorted.read((char*)node_read_buffer.data(), node_read_buffer.size() * sizeof(Node));

        cout << "Read " << node_order_sorted.gcount() << " bytes (" << (node_order_sorted.gcount() / sizeof(Node)) << " elements)" << endl;
        node_read_buffer.resize(node_order_sorted.gcount() / sizeof(Node));
        for (auto node : node_read_buffer) {
            // cout << node.height << endl;
            if (!flow_lengths.inGrid(node.x + dx[node.dir], node.y + dy[node.dir])) {
                // SKipping flow out of grid
                continue;
            }
            if (dem(node.x, node.y) <= dem(node.x + dx[node.dir], node.y + dy[node.dir])) {
                cout << "Invalid sort " << dem(node.x, node.y) << " " << dem(node.x + dx[node.dir], node.y + dy[node.dir]) << endl;
            }
            assert(dem(node.x, node.y) >= dem(node.x + dx[node.dir], node.y + dy[node.dir]));
            auto current = flow_lengths(node.x, node.y);
            auto& next = flow_lengths(node.x + dx[node.dir], node.y + dy[node.dir]);
            next = max(next, current + 1);
            record = max(next, record);
        }
    } while(node_order_sorted);
    node_order_sorted.close();
    */

    cout << "Record length is " << record << endl;
    return flow_lengths;
}

int main (int argc, char** argv) {
    if (argc != 4) {
        cerr << "Usage: ./flow_length layoutfile_dem layoutfile_directions outpath" << endl;
        cerr << "outpath must include a %f" << endl;
        return 1;
    }

    assert(string(argv[3]).find_last_of("/") != string::npos);
    assert(string(argv[3]).find("%f") != string::npos);

    int cache_size = 16000*6;
    // cout << "Reading tile metadata..." << flush;
    // A2Array2D<float> dem(argv[1], cache_size);
    // cout << " done" << endl;

    cout << "Reading tile metadata..." << flush;
    A2Array2D<flowdir_t> directions(argv[2], cache_size);
    cout << " done" << endl;

    A2Array2D<int> output = calculate_flow(directions, cache_size);

    auto foldername = string(argv[3]).substr(0, string(argv[3]).find_last_of("/")).c_str();
    mkdir(foldername, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    cout << "Saving final output" << endl;
    output.saveGDAL(argv[3]);
    cout << "Saving unified gdal" << endl;
    output.saveUnifiedGDAL("flow_length.gdal");
    
    return 0;
}