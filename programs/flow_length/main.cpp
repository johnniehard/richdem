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


//Valid flowdirs are in the range 0-8, inclusive. 0 is the center and 1-8,
//inclusive, are the 8 cells around the central cell. An extra value is needed
//to indicate NoData. Therefore, uint8_t is appropriate.
typedef uint8_t flowdir_t;

A2Array2D<float> calculate_flow(A2Array2D<flowdir_t> &directions, int cache_size) {
    A2Array2D<int> in_count("tmp/in_count/", directions.stdTileWidth(), directions.stdTileHeight(), directions.widthInTiles(), directions.heightInTiles(), cache_size);
    in_count.setAll(0);

    cout << "Counting arrows pointing in" << endl;
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

    cout << endl;
    cout << "Checking for leaf nodes" << endl;
    queue<tuple<int,int>> que;
    for (int y = 0; y < directions.height(); y++) {
        cout << "\r" << (int)(100*(float)(y+1)/directions.height()) << "%" << flush;

        for (int x = 0; x < directions.width(); x++) {
            if (in_count(x, y) == 0) que.push({x, y});
        }
    }

    cout << endl;
    cout << "Calculating flow lengths" << endl;
    A2Array2D<float> flow_lengths("tmp/flow_length/", directions.stdTileWidth(), directions.stdTileHeight(), directions.widthInTiles(), directions.heightInTiles(), cache_size);
    flow_lengths.copy_metadata_from(directions);
    flow_lengths.setAll(0);

    float pixelSize = directions.getCellLength();

    float diagLength = sqrt(2.0f) * pixelSize;
    float straightLength = pixelSize;

    float record = 0;
    int64_t done = 0;
    while(!que.empty()) {
        done++;
        if ((done % (1024*128)) == 0){
            cout << "\r" << (int)(( done / (double)((int64_t)directions.width()*directions.height())) * 100) << "%. Loops done:" << done << flush;
        }
        
        int x, y;
        tie(x, y) = que.front();
        que.pop();
        auto dir = directions(x, y);
        auto nx = x + dx[dir];
        auto ny = y + dy[dir];
        if (!flow_lengths.in_grid(nx, ny)) continue;

        auto current = flow_lengths(x, y);
        auto& next = flow_lengths(nx, ny);
        float step_distance = n_diag[dir] ? diagLength : straightLength;
        next = max(next, current + step_distance);
        record = max(next, record);
        auto& next_in = in_count(nx, ny);
        next_in--;
        if (next_in == 0) que.push({nx, ny});
    }

    cout << endl;
    cout << "Record length is " << record << endl;
    return flow_lengths;
}

int main (int argc, char** argv) {
    if (argc != 3) {
        cerr << "Usage: ./flow_length layoutfile_directions outpath" << endl;
        cerr << "outpath must include a %f" << endl;
        return 1;
    }

    assert(string(argv[2]).find_last_of("/") != string::npos);
    assert(string(argv[2]).find("%f") != string::npos);

    int cache_size = 16000*6;

    cout << "Reading tile metadata..." << flush;
    A2Array2D<flowdir_t> directions(argv[1], cache_size);
    cout << " done" << endl;

    auto output = calculate_flow(directions, cache_size);

    auto foldername = string(argv[2]).substr(0, string(argv[2]).find_last_of("/")).c_str();
    mkdir(foldername, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    cout << "Saving final output" << endl;
    output.saveGDAL(argv[2]);

    cout << "Saving unified gdal (only for debug)" << endl;
    output.saveUnifiedGDAL("flow_length.gdal");
    
    return 0;
}