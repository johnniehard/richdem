#include "richdem/tiled/A2Array2D.hpp"
#include <iostream>
#include <vector>
#include <chrono>
#include <thread>
#include <string>
#include <queue>
#include <channel.h>
#include <mutex>
#include <atomic>
#include <future>
#include <sys/stat.h>

using namespace std;
using namespace richdem;

enum class TileState {
    Unloaded,
    Loading,
    Loaded,
    Saving,
};

template<class T>
struct Tile {
    vector<T> elevation;
    vector<T> breachingElevation;
    vector<int> distanceFromSink;
    vector<uint8_t> parent;

    string filename = "";
    int x = -1, y = -1;
    int tileSize = -1;
    atomic<TileState> state;
    atomic<bool> changed;
    mutex tileMutex;
    T noData;
    atomic<int> referenceCount;

    Tile(int x, int y, int tileSize, string filename) : filename(filename), x(x), y(y), tileSize(tileSize), state(TileState::Unloaded), changed(false), referenceCount(0) {}

    void loadFrom(A2Array2D<T>& grid, int xOffset, int yOffset) {
        lock_guard<mutex> guard(tileMutex);
        assert(state == TileState::Unloaded);

        elevation = vector<T>(tileSize * tileSize);
        breachingElevation = vector<T>(tileSize * tileSize, numeric_limits<T>::max());
        distanceFromSink = vector<int>(tileSize * tileSize, numeric_limits<int>::max());
        parent = vector<uint8_t>(tileSize * tileSize, numeric_limits<uint8_t>::max());
        noData = grid.noData();

        for (int y = 0; y < tileSize; y++) {
            for (int x = 0; x < tileSize; x++) {
                if (grid.inGrid(xOffset + x, yOffset + y)) {
                    elevation[y * tileSize + x] = grid(xOffset + x, yOffset + y);
                } else {
                    elevation[y * tileSize + x] = noData;
                }
            }
        }

        state = TileState::Loaded;
    }

    void load() {
        if (state == TileState::Loaded) return;
        lock_guard<mutex> guard(tileMutex);

        assert(state == TileState::Unloaded);

        state = TileState::Loading;

        std::ifstream fin(filename, std::ios::in | std::ios::binary);
        assert(fin.good());

        elevation.resize(tileSize*tileSize);
        breachingElevation.resize(tileSize*tileSize);
        distanceFromSink.resize(tileSize*tileSize);
        parent.resize(tileSize*tileSize);

        assert((int)parent.size() == tileSize*tileSize);

        fin.read(reinterpret_cast<char*>(elevation.data()), tileSize*tileSize*sizeof(T));
        fin.read(reinterpret_cast<char*>(breachingElevation.data()), tileSize*tileSize*sizeof(T));
        fin.read(reinterpret_cast<char*>(distanceFromSink.data()), tileSize*tileSize*sizeof(int));
        fin.read(reinterpret_cast<char*>(parent.data()), tileSize*tileSize*sizeof(uint8_t));

        state = TileState::Loaded;
    }

    void save() {
        lock_guard<mutex> guard(tileMutex);
        assert(state == TileState::Loaded);
        state = TileState::Saving;

        std::ofstream fout(filename, std::ios::out | std::ios::binary);
        assert(fout.good());
        assert((int)elevation.size() == tileSize*tileSize);

        fout.write(reinterpret_cast<char*>(elevation.data()), tileSize*tileSize*sizeof(T));
        fout.write(reinterpret_cast<char*>(breachingElevation.data()), tileSize*tileSize*sizeof(T));
        fout.write(reinterpret_cast<char*>(distanceFromSink.data()), tileSize*tileSize*sizeof(int));
        fout.write(reinterpret_cast<char*>(parent.data()), tileSize*tileSize*sizeof(uint8_t));
        elevation = vector<T>();
        breachingElevation = vector<T>();
        distanceFromSink = vector<int>();
        parent = vector<uint8_t>();

        state = TileState::Unloaded;
    }
};

template<class T>
struct QuadGrid {
    vector<shared_ptr<Tile<T>>> data;
    int tileSize;

    QuadGrid(int tileSize, vector<shared_ptr<Tile<T>>> data) : data(data), tileSize(tileSize) {}

    int width() const {
        return tileSize*2;
    }

    int height() const {
        return tileSize*2;
    }

    bool contains(int x, int y) const {
        return x >= 0 && y >= 0 && x < tileSize*2 && y < tileSize*2;
    }

    bool isEdgeCell(int x, int y) const {
        int x0 = data[0]->x;
        int y0 = data[0]->y;
        if (x == 0 && x0 == 0) return true;
        if (y == 0 && y0 == 0) return true;
        // TODO: Need to identify max-x and max-y edges
        return false;
    }

    T& elevation(int x, int y) {
        int x0 = x;
        int y0 = y;
        int tileIndex = 0;
        if (y >= tileSize) {
            y0 -= tileSize;
            tileIndex += 2;
        }
        if (x >= tileSize) {
            x0 -= tileSize;
            tileIndex += 1;
        }
        return data[tileIndex]->elevation[y0 * tileSize + x0];
    }

    T elevation(int x, int y) const {
        T& v = (*this)(x, y);
        return v;
    }

    int& distanceFromSink(int x, int y) {
        int x0 = x;
        int y0 = y;
        int tileIndex = 0;
        if (y >= tileSize) {
            y0 -= tileSize;
            tileIndex += 2;
        }
        if (x >= tileSize) {
            x0 -= tileSize;
            tileIndex += 1;
        }
        return data[tileIndex]->distanceFromSink[y0 * tileSize + x0];
    }

    T& breachingElevation(int x, int y) {
        int x0 = x;
        int y0 = y;
        int tileIndex = 0;
        if (y >= tileSize) {
            y0 -= tileSize;
            tileIndex += 2;
        }
        if (x >= tileSize) {
            x0 -= tileSize;
            tileIndex += 1;
        }
        return data[tileIndex]->breachingElevation[y0 * tileSize + x0];
    }

    uint8_t& parent(int x, int y) {
        int x0 = x;
        int y0 = y;
        int tileIndex = 0;
        if (y >= tileSize) {
            y0 -= tileSize;
            tileIndex += 2;
        }
        if (x >= tileSize) {
            x0 -= tileSize;
            tileIndex += 1;
        }
        return data[tileIndex]->parent[y0 * tileSize + x0];
    }
};

template<class T>
struct ElevationPixel {
    int x, y;
    T elevation;
    int distanceFromSink;

    ElevationPixel (int x, int y, T elevation, int distanceFromSink) : x(x), y(y), elevation(elevation), distanceFromSink(distanceFromSink) {}

    bool operator< (const ElevationPixel& other) const {
        return elevation < other.elevation || (elevation == other.elevation && distanceFromSink < other.distanceFromSink);
    }

    bool operator> (const ElevationPixel& other) const {
        return elevation > other.elevation || (elevation == other.elevation && distanceFromSink > other.distanceFromSink);
    }
};

template<class T>
void loadThread(vector<vector<shared_ptr<Tile<T>>>> grid, channel<QuadGrid<T>>& outChannel, channel<shared_ptr<Tile<T>>>& saveChannel) {
    int gridWidth = grid[0].size();
    int gridHeight = grid.size();
    int tileSize = grid[0][0]->tileSize;

    vector<QuadGrid<T>> loadOrder;
    // Down-Right
    for (int x = 0; x < gridWidth - 1; x++) {
        for (int y = 0; y < gridHeight - 1; y++) {
            loadOrder.push_back(QuadGrid<T>(tileSize, { grid[y][x], grid[y][x+1], grid[y+1][x], grid[y+1][x+1] }));
        }
    }

    QuadGrid<T> lastQuad = QuadGrid<T>(0, {});
    for (auto& quad : loadOrder) {
        for (auto tile : quad.data) {
            tile->load();
            tile->referenceCount++;
        }
        for (auto tile : lastQuad.data) {
            if (--tile->referenceCount == 0) {
                saveChannel.put(tile);
            }
        }

        outChannel.put(quad);
        lastQuad = quad;
    }

    outChannel.close();
}

template<class T>
void runDijkstra(QuadGrid<T>& quad) {
    priority_queue<ElevationPixel<T>, vector<ElevationPixel<T>>, std::greater<ElevationPixel<T>>> que;
    vector<vector<uint8_t>> visited (quad.height(), vector<uint8_t>(quad.width()));
    auto sinkValue = quad.data[0]->noData;

    int hash0 = 0;
    for (int y = 0; y < quad.height(); y++) {
        for (int x = 0; x < quad.width(); x++) {
            hash0 = 31 * hash0 ^ quad.parent(x, y);
            hash0 = 31 * hash0 ^ quad.distanceFromSink(x, y);
            hash0 = 31 * hash0 ^ (int)(1000*quad.breachingElevation(x, y));
        }
    }
    cout << "Hash " << hash0 << endl;

    // TODO: Add when tile is created
    int sinks = 0;
    for (int y = 0; y < quad.height(); y++) {
        for (int x = 0; x < quad.width(); x++) {
            int distanceFromSink = quad.distanceFromSink(x, y);
            if (distanceFromSink == numeric_limits<int>::max() && quad.elevation(x, y) == sinkValue) {
                quad.breachingElevation(x, y) = 0;
                quad.distanceFromSink(x, y) = 0;
                sinks++;
            }

            if (quad.isEdgeCell(x, y)) {
                quad.breachingElevation(x, y) = quad.elevation(x, y);
                quad.distanceFromSink(x, y) = 0;
                sinks++;
            }

            if (quad.distanceFromSink(x, y) != numeric_limits<int>::max()) {
                //que.emplace(x, y, quad.elevation(x, y), quad.distanceFromSink(x, y));
                visited[y][x] = 1;
            }
        }
    }
    

    cout << "Added " << sinks << " sinks" << endl;

    int c = 0;
    for (int y = 0; y < quad.height(); y++) {
        for (int x = 0; x < quad.width(); x++) {
            if (visited[y][x] == 1) {
                int distanceFromSink = quad.distanceFromSink(x, y);
                T breachElevation = quad.breachingElevation(x, y);
                c++;
                
                // This pixel has been visited
                for (int n = 0; n < 8; n++) {
                    const int nx = x + dx[n+1];
                    const int ny = y + dy[n+1];

                    if (quad.contains(nx, ny)) {
                        auto nelev = quad.elevation(nx, ny);
                        auto potentialBreach = max(breachElevation, nelev);
                        auto& nBreachingElevation = quad.breachingElevation(nx, ny);
                        if (potentialBreach < nBreachingElevation || (potentialBreach == nBreachingElevation && distanceFromSink + 1 < quad.distanceFromSink(nx, ny))) {
                            nBreachingElevation = potentialBreach;
                            quad.parent(nx, ny) = n;
                            quad.distanceFromSink(nx, ny) = distanceFromSink + 1;
                            que.emplace(nx, ny, nelev, distanceFromSink + 1);
                            visited[y][x] = 0;
                        }
                    }
                }
            }
        }
    }

    cout << "Added " << que.size() << " seeds" << endl;

    // Mark the tiles as being changed (for debugging purposes) if we did change anything

    bool changed = false;

    int k = 0;
    while(!que.empty()) {
        auto pixel = que.top();
        que.pop();
        k++;

        if ((k % 1000000) == 0) cout << k << " " << que.size() << " " << pixel.elevation << endl;

        uint8_t isVisited = visited[pixel.y][pixel.x];

        // Ignore if we have already visited this pixel
        if (isVisited == 2) continue;
        visited[pixel.y][pixel.x] = 2;
        auto breachElevation = quad.breachingElevation(pixel.x, pixel.y);
        int distanceFromSink = quad.distanceFromSink(pixel.x, pixel.y);

        for (int n = 0; n < 8; n++) {
            const int nx = pixel.x + dx[n+1];
            const int ny = pixel.y + dy[n+1];

            if (quad.contains(nx, ny)) {
                auto nelev = quad.elevation(nx, ny);
                auto potentialBreach = max(breachElevation, nelev);
                auto& nBreachingElevation = quad.breachingElevation(nx, ny);
                if (potentialBreach < nBreachingElevation || (potentialBreach == nBreachingElevation && distanceFromSink + 1 < quad.distanceFromSink(nx, ny))) {
                    assert(visited[ny][nx] != 2);
                    nBreachingElevation = potentialBreach;
                    quad.parent(nx, ny) = n;
                    quad.distanceFromSink(nx, ny) = distanceFromSink + 1;
                    que.emplace(nx, ny, nelev, distanceFromSink + 1);
                    changed = true;
                }
            }
        }
    }

    int hash = 0;
    for (int y = 0; y < quad.height(); y++) {
        for (int x = 0; x < quad.width(); x++) {
            hash = 31 * hash ^ quad.parent(x, y);
            hash = 31 * hash ^ quad.distanceFromSink(x, y);
            hash = 31 * hash ^ (int)(1000*quad.breachingElevation(x, y));
        }
    }
    cout << "Hash " << hash << endl;

    for (auto tile : quad.data) tile->changed = tile->changed | changed;
}

template<class T>
void dijkstraThread(channel<QuadGrid<T>>& workChannel, channel<shared_ptr<Tile<T>>>& saveChannel) {
    QuadGrid<T> quad = QuadGrid<T>(0, {});
    while(workChannel.get(quad)) {
        runDijkstra(quad);

        for (auto tile : quad.data) {
            if (--tile->referenceCount == 0) {
                saveChannel.put(tile);
            }
        }
    }

    saveChannel.close();
}

template<class T>
void saveThread (channel<shared_ptr<Tile<T>>>& saveChannel) {
    shared_ptr<Tile<T>> tile;
    while(saveChannel.get(tile)) {
        tile->save();
    }
}

template<class T>
vector<vector<shared_ptr<Tile<T>>>> createTiles(string path, A2Array2D<T>& grid, int tileSize, channel<shared_ptr<Tile<T>>>& saveChannel) {

    //assert((grid.width() % tileSize) == 0);
    //assert((grid.height() % tileSize) == 0);

    int gridWidth = (grid.width()+tileSize-1) / tileSize;
    int gridHeight = (grid.height()+tileSize-1) / tileSize;
    vector<vector<shared_ptr<Tile<T>>>> result(gridHeight, vector<shared_ptr<Tile<T>>>(gridWidth));

    for (int y = 0; y < gridHeight; y++) {
        for (int x = 0; x < gridWidth; x++) {
            auto tile = make_shared<Tile<T>>(x, y, tileSize, "tmp/tiles/" + to_string(x) + "_" + to_string(y) + ".bin");
            result[y][x] = tile;
            tile->loadFrom(grid, x * gridWidth, y * gridWidth);
            saveChannel.put(tile);
        }
    }

    saveChannel.close();
    return result;
}

void testProducer(channel<float>& c) {
    for (int i =0; i < 1000; i++) {
        cout << "Putting " << i << endl;
        c.put(i);
    }
    c.close();
}

void testConsumer(channel<float>& c) {
    for (int i =0; i < 1000; i++) {
        float v;
        if (c.get(v)) {
            cout << "Got " << v << endl;
        } else {
            assert(false);
        }
        //this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    float v2;
    if (c.get(v2)) {
        cout << "ERROR" << endl;
    }
}

void minitor(std::function<bool()> until, vector<vector<shared_ptr<Tile<float>>>>& grid) {
    while(!until()) {
        // Move cursor to top left of screen
        cout << "\x1b[1;1H";
        for (int y = 0; y < (int)grid.size(); y++) {
            for (int x = 0; x < (int)grid[y].size(); x++) {
                auto tile = grid[y][x];
                string c = "_";
                string color = "";
                switch(tile->state) {
                case TileState::Unloaded:
                    c = "#";
                    color = "50;50;50";
                    break;
                case TileState::Loading:
                    c = "#";
                    color = "55;126;184";
                    break;
                case TileState::Loaded:
                    c = "#";
                    color = "77;175;74";
                    break;
                case TileState::Saving:
                    c = "#";
                    color = "228;26;28";
                    break;
                }
                cout << "\x1b[38;2;" << color << "m" << c << "\x1b[0m";
            }
            // Clear from cursor to end of line
            cout << "\x1b[0K";
            cout << "\n";
        }
        // Clear from the cursor to end of screen
        cout << "\x1b[0J";
        cout << std::flush;
        this_thread::sleep_for(std::chrono::milliseconds(16));
    }
}

int main () {
    // channel<float> c(8);
    // thread a(testProducer, std::ref(c));
    // thread b(testConsumer, std::ref(c));

    // a.join();
    // b.join();
    // return 0;

    int tileSize = 4096;

    mkdir("tmp", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir("tmp/tiles", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


    vector<vector<shared_ptr<Tile<float>>>> grid;
    {
        A2Array2D<float> gdalGrid("out_tiler/layout.layout", 128);
        channel<shared_ptr<Tile<float>>> saveChannel(16);
        auto future = std::async(createTiles<float>, "", std::ref(gdalGrid), tileSize, std::ref(saveChannel));
        thread saveT(saveThread<float>, std::ref(saveChannel));

        
        saveT.join();

        grid = future.get();
    }

    while(true) {
        channel<QuadGrid<float>> workChannel(8);
        channel<shared_ptr<Tile<float>>> saveChannel(8);

        thread loadT (loadThread<float>, grid, std::ref(workChannel), std::ref(saveChannel));
        thread dijkstraT (dijkstraThread<float>, std::ref(workChannel), std::ref(saveChannel));
        thread saveT(saveThread<float>, std::ref(saveChannel));
        
        //minitor([&]() { return saveChannel.is_closed() && saveChannel.size() == 0; }, grid);
        loadT.join();
        dijkstraT.join();
        saveT.join();

        bool changed = false;
        for (int y = 0; y < (int)grid.size(); y++) {
            for (int x = 0; x < (int)grid[y].size(); x++) {
                auto tile = grid[y][x];
                changed |= tile->changed;
                tile->changed = false;
            }
        }

        if (!changed) break;
    }

    // thread load(loadThread, grid, 
    cout << "Completed" << endl;
    return 0;
}