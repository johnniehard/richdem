/**
  @file
  @brief Defines structures for addressing grid cells and associated queues.

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_grid_hpp_
#define _richdem_grid_hpp_

#include <vector>
#include <queue>
#include <cmath>
#include <functional>

namespace richdem {

/// Stores the (x,y) coordinates of a grid cell
class GridCell {
 public:
  int x; ///< Grid cell's x-coordinate
  int y; ///< Grid cell's y-coordinate
  /// Initiate the grid cell without coordinates; should generally be avoided
  GridCell(){}
  /// Initiate the grid cell to the coordinates (x0,y0)
  GridCell(int x, int y):x(x),y(y){}
};


///@brief Stores the (x,y,z) coordinates of a grid cell; useful for priority sorting with \ref GridCellZ_pq.
template<class elev_t>
class GridCellZ : public GridCell {
 public:
  elev_t z;         ///< Grid cell's z-coordinate
  GridCellZ(int x, int y, elev_t z): GridCell(x,y), z(z) {}
  GridCellZ(){}
  //TODO: Add tests of isnan
  bool isnan() const { return false; }
  bool operator> (const GridCellZ<elev_t>& a) const { return z>a.z; }
};

///@brief An (x,y,z) cell with NaNs taken as infinitely small numbers.
template<>
class GridCellZ<double> : public GridCell {
 public:
  double z;         ///< Grid cell's z-coordinate
  GridCellZ(int x, int y, double z): GridCell(x,y), z(z) {}
  GridCellZ(){}
  //TODO: Add tests of isnan
  bool isnan() const { return std::isnan(z); }
  bool operator< (const GridCellZ<double>& a) const { return ( std::isnan(z) && !std::isnan(a.z)) || z< a.z; }
  bool operator> (const GridCellZ<double>& a) const { return (!std::isnan(z) &&  std::isnan(a.z)) || z> a.z; }
  bool operator>=(const GridCellZ<double>& a) const { return ( std::isnan(z) &&  std::isnan(a.z)) || (!std::isnan(z) &&  std::isnan(a.z)) || z>=a.z; }
  bool operator<=(const GridCellZ<double>& a) const { return ( std::isnan(z) &&  std::isnan(a.z)) || ( std::isnan(z) && !std::isnan(a.z)) || z<=a.z; }
  bool operator==(const GridCellZ<double>& a) const { return ( std::isnan(z) &&  std::isnan(a.z)) || z==a.z; }
  bool operator!=(const GridCellZ<double>& a) const { return !(std::isnan(z) &&  std::isnan(a.z)) && z!=a.z; }
};

///@brief An (x,y,z) cell with NaNs taken as infinitely small numbers.
template<>
class GridCellZ<float>: public GridCell {
 public:
  float z;         ///< Grid cell's z-coordinate
  GridCellZ(int x, int y, float z): GridCell(x,y), z(z) {}
  GridCellZ(){}
  //TODO: Add tests of isnan
  ///Compare cells based on elevation. (TODO: Distribute)
  bool isnan() const { return std::isnan(z); }
  bool operator< (const GridCellZ<float>& a) const { return ( std::isnan(z) && !std::isnan(a.z)) || z< a.z; }
  bool operator> (const GridCellZ<float>& a) const { return (!std::isnan(z) &&  std::isnan(a.z)) || z> a.z; }
  bool operator>=(const GridCellZ<float>& a) const { return ( std::isnan(z) &&  std::isnan(a.z)) || (!std::isnan(z) &&  std::isnan(a.z)) || z>=a.z; }
  bool operator<=(const GridCellZ<float>& a) const { return ( std::isnan(z) &&  std::isnan(a.z)) || ( std::isnan(z) && !std::isnan(a.z)) || z<=a.z; }
  bool operator==(const GridCellZ<float>& a) const { return ( std::isnan(z) &&  std::isnan(a.z)) || z==a.z; }
  bool operator!=(const GridCellZ<float>& a) const { return !(std::isnan(z) &&  std::isnan(a.z)) && z!=a.z; }
};


///@brief Stores the (x,y,z) coordinates of a grid cell and a priority indicator k; used by \ref GridCellZk_pq.
template<class elev_t>
class GridCellZk : public GridCellZ<elev_t> {
  public:
    int k;           ///< Used to store an integer to make sorting stable
    GridCellZk(int x, int y, elev_t z, int k): GridCellZ<elev_t>(x,y,z), k(k) {}
    GridCellZk(){}
    //TODO: Is it possible to do this relying on inheriting the std::isnan checks from the GridCellZ specialization?
    // bool operator< (const GridCellZk<elev_t>& a) const { return GridCellZ<elev_t>::z< a.z || ( GridCellZ<elev_t>::isnan() && !a.isnan()) || (GridCellZ<elev_t>::z==a.z && k<a.k) || (GridCellZ<elev_t>::isnan() && a.isnan() && k<a.k); }
    // bool operator> (const GridCellZk<elev_t>& a) const { return GridCellZ<elev_t>::z> a.z || (!GridCellZ<elev_t>::isnan() &&  a.isnan()) || (GridCellZ<elev_t>::z==a.z && k>a.k) || (GridCellZ<elev_t>::isnan() && a.isnan() && k>a.k); }

    bool operator< (const GridCellZk<elev_t>& a) const { return GridCellZ<elev_t>::z < a.z || (GridCellZ<elev_t>::z==a.z && k<a.k); }
    bool operator> (const GridCellZk<elev_t>& a) const { return GridCellZ<elev_t>::z > a.z || (GridCellZ<elev_t>::z == a.z && k>a.k); }
};





///@brief A priority queue of GridCellZ, sorted by ascending height
template<typename elev_t>
using GridCellZ_pq = std::priority_queue<GridCellZ<elev_t>, std::vector<GridCellZ<elev_t> >, std::greater<GridCellZ<elev_t> > >;


///@brief A priority queue of GridCellZk, sorted by ascending height or, if heights are equal, by the order of insertion.
template<typename T>
class GridCellZk_pq : public std::priority_queue<GridCellZk<T>, std::vector< GridCellZk<T> >, std::greater<GridCellZk<T> > > {
 private:
  uint64_t count = 0;
 public:
  void push(){ //TODO: Is there a way to stop compilation, but only if this function is used
    throw std::runtime_error("push() to GridCellZk_pq is not allowed!");
  }
  void emplace(int x, int y, T z){
    std::priority_queue<GridCellZk<T>, std::vector< GridCellZk<T> >, std::greater<GridCellZk<T> > >::emplace(x,y,z,++count);
  }
};

///@brief Stores the (x,y,z) coordinates of a grid cell and a priority indicator k; used by \ref GridCellZk_pq.
template<class elev_t>
class GridCellZkB {
  public:
    elev_t z;
    int x, y;
    uint32_t k;           ///< Used to store an integer to make sorting stable
    int backlink_dir;

    GridCellZkB(int x, int y, elev_t z, uint32_t k, int backlink_dir): z(z), x(x), y(y), k(k), backlink_dir(backlink_dir) {}
    GridCellZkB(){}
    //TODO: Is it possible to do this relying on inheriting the std::isnan checks from the GridCellZ specialization?
    // bool operator< (const GridCellZk<elev_t>& a) const { return GridCellZ<elev_t>::z< a.z || ( GridCellZ<elev_t>::isnan() && !a.isnan()) || (GridCellZ<elev_t>::z==a.z && k<a.k) || (GridCellZ<elev_t>::isnan() && a.isnan() && k<a.k); }
    // bool operator> (const GridCellZk<elev_t>& a) const { return GridCellZ<elev_t>::z> a.z || (!GridCellZ<elev_t>::isnan() &&  a.isnan()) || (GridCellZ<elev_t>::z==a.z && k>a.k) || (GridCellZ<elev_t>::isnan() && a.isnan() && k>a.k); }

    bool operator< (const GridCellZkB<elev_t>& a) const { return z < a.z || (z==a.z && k<a.k); }
    bool operator> (const GridCellZkB<elev_t>& a) const { return z > a.z || (z == a.z && k>a.k); }
};

///@brief A priority queue of GridCellZkB, sorted by ascending height or, if heights are equal, by the order of insertion.
template<typename T>
class GridCellZkB_pq : public std::priority_queue<GridCellZkB<T>, std::vector< GridCellZkB<T> >, std::greater<GridCellZkB<T> > > {
 private:
  uint32_t count = 0;
 public:
  void push(){ //TODO: Is there a way to stop compilation, but only if this function is used
    throw std::runtime_error("push() to GridCellZk_pq is not allowed!");
  }
  void emplace(int x, int y, T z, int backlink_dir){
    std::priority_queue<GridCellZkB<T>, std::vector< GridCellZkB<T> >, std::greater<GridCellZkB<T> > >::emplace(x, y, z, ++count, backlink_dir);
  }
};

}

#endif
