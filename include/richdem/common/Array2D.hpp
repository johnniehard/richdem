#ifndef _array_2d_hpp_
#define _array_2d_hpp_

#include "gdal_priv.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <typeinfo>
#include <stdexcept>
#include <unordered_set> //For printStamp
#include "richdem/common/constants.hpp"

//These enable compression in the loadNative() and saveNative() methods
#ifdef WITH_COMPRESSION
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#endif

/**
  @brief  Determine data type of a GDAL file's first layer
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]  filename   Filename of file whose type should be determined
*/
GDALDataType peekGDALType(const std::string &filename) {
  GDALAllRegister();
  GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
  assert(fin!=NULL);

  GDALRasterBand *band   = fin->GetRasterBand(1);
  GDALDataType data_type = band->GetRasterDataType();

  GDALClose(fin);

  return data_type;
}


/**
  @brief  Retrieve height, width, NoData, and geotransform from a GDAL file
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]  filename   GDAL file to peek at
  @param[out] height     Height of the raster in cells
  @param[out] width      Width of the raster in cells
  @param[out] no_data    Value of the raster's no_data property
  @param[out] geo_trans  Returns the SIX elements of the raster's geotransform
*/
template<class T>
void getGDALHeader(
  const   std::string &filename,
  int32_t &height,
  int32_t &width,
  T       &no_data,
  double  geotransform[6]
){
  GDALAllRegister();
  GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
  assert(fin!=NULL);

  GDALRasterBand *band   = fin->GetRasterBand(1);

  height  = band->GetYSize();
  no_data = band->GetNoDataValue();
  width   = band->GetXSize();

  fin->GetGeoTransform(geotransform);

  GDALClose(fin);
}


/**
  @brief  Retrieve height, width, data type, and geotransform from a GDAL file
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]  filename   GDAL file to peek at
  @param[out] height     Height of the raster in cells
  @param[out] width      Width of the raster in cells
  @param[out] dtype      Data type of the file in question
  @param[out] geo_trans  Returns the SIX elements of the raster's geotransform
*/
void getGDALDimensions(
  const   std::string &filename,
  int32_t &height,
  int32_t &width,
  GDALDataType &dtype,
  double geotransform[6]
){
  GDALAllRegister();
  GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
  assert(fin!=NULL);

  GDALRasterBand *band = fin->GetRasterBand(1);
  
  dtype = band->GetRasterDataType();

  if(geotransform!=NULL && fin->GetGeoTransform(geotransform)!=CE_None){
    std::cerr<<"Error getting geotransform from '"<<filename<<"'!"<<std::endl;
    throw std::runtime_error("Could not get geotransform!");
  }

  height  = band->GetYSize();
  width   = band->GetXSize();

  GDALClose(fin);
}


/**
  @brief  Convert Array2D or any other template to its GDAL data type
  @author Richard Barnes (rbarnes@umn.edu)

  @return The GDAL datatype of T
*/
template<class T>
GDALDataType NativeTypeToGDAL() {
  if(typeid(T)==typeid(uint8_t))
    return GDT_Byte;
  else if(typeid(T)==typeid(uint16_t))
    return GDT_UInt16;
  else if(typeid(T)==typeid(int16_t))
    return GDT_Int16;
  else if(typeid(T)==typeid(uint32_t))
    return GDT_UInt32;
  else if(typeid(T)==typeid(int32_t))
    return GDT_Int32;
  else if(typeid(T)==typeid(float))
    return GDT_Float32;
  else if(typeid(T)==typeid(double))
    return GDT_Float64;
  else {
    std::cerr<<"Could not map native type '"<<typeid(T).name()<<"' to GDAL type! (Use `c++filt -t` to decode.)"<<std::endl;
    throw std::runtime_error("Could not map native data type to GDAL type!");
  }
  return GDT_Unknown;
}




/**
  @brief  Class to hold and manipulate GDAL and native rasters
  @author Richard Barnes (rbarnes@umn.edu)

  Array2D manages a two-dimensional raster dataset. Passed a request to load
  such data, it peeks at the file header and can either load data on
  construction or wait until a later point. It can also offload data to disk.

  Array2D permits simple copy construction as well as templated copies, which
  transfer projections and geotransforms, but not the actual data. This is
  useful for say, create a flow directions raster which is homologous to a DEM.
*/
template<class T>
class Array2D {
 public:
  std::string filename;             ///TODO
  std::string basename;             ///Filename without path or extension
  std::vector<double> geotransform; ///Geotransform of the raster
  std::string projection;           ///Projection of the raster

 private:
  template<typename> friend class Array2D;

  std::vector<T> data;              ///Holds the raster data in a 1D array
                                    ///this improves caching versus a 2D array

  T   no_data;                      ///NoData value of the raster
  int64_t num_data_cells = -1;      ///Number of cells which are not NoData

  int32_t view_width;               ///Height of raster in cells
  int32_t view_height;              ///Width of raster in cells

  ///@{ A rectangular subregion of a larger raster can be extracted. These
  ///   variables store the offsets of this subregion in case the subregion
  ///   needs to be saved into a raster with other subregions
  int32_t view_xoff;
  int32_t view_yoff;
  ///@}
  
  bool from_cache;  ///If TRUE, loadData() loads data from the cache assuming 
                    ///the Native format. Otherwise, it assumes it is loading
                    ///from a GDAL file.

  ///TODO
  void loadGDAL(const std::string &filename, int32_t xOffset=0, int32_t yOffset=0, int32_t part_width=0, int32_t part_height=0, bool exact=false, bool load_data=true){
    assert(empty());

    from_cache = false;

    this->filename = filename;

    GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
    if(fin==NULL){
      std::cerr<<"Could not open file '"<<filename<<"'!"<<std::endl;
      throw std::runtime_error("Could not open a GDAL file!");
    }

    geotransform.resize(6);
    if(fin->GetGeoTransform(geotransform.data())!=CE_None){
      std::cerr<<"Error getting geotransform from '"<<filename<<"'!"<<std::endl;
      throw std::runtime_error("Error getting geotransform!");
    }

    const char* projection_string=fin->GetProjectionRef();
    projection = std::string(projection_string);

    GDALRasterBand *band = fin->GetRasterBand(1);

    int32_t total_width  = band->GetXSize();         //Returns an int
    int32_t total_height = band->GetYSize();         //Returns an int
    no_data              = band->GetNoDataValue();

    if(exact && (total_width-xOffset!=part_width || total_height-yOffset!=part_height))
      throw std::runtime_error("Tile dimensions did not match expectations!");

    //TODO: What's going on here?

    if(xOffset+part_width>=total_width)
      part_width  = total_width-xOffset;
    if(yOffset+part_height>=total_height)
      part_height = total_height-yOffset;

    if(part_width==0)
      part_width = total_width;
    view_width = part_width;

    if(part_height==0)
      part_height = total_height;
    view_height = part_height;

    view_xoff = xOffset;
    view_yoff = yOffset;

    GDALClose(fin);

    //std::cerr<<"Allocating: "<<view_height<<" rows by "<<view_width<<" columns"<<std::endl;
    if(load_data)
      loadData();
  }

  ///Returns the GDAL data type of the Array2D template type
  GDALDataType myGDALType() const {
    return NativeTypeToGDAL<T>();
  }

  /**
    @brief Saves raster to a simply-structure file on disk, possibly using
           compression.

    @post  Using loadData() after running this function will result in data
           being loaded from the cache, rather than the original file (if any).
  */
  void saveToCache(const std::string &filename){
    std::fstream fout;

    from_cache     = true;
    this->filename = filename;

    fout.open(filename, std::ios_base::binary | std::ios_base::out | std::ios::trunc);
    if(!fout.good()){
      std::cerr<<"Failed to open file '"<<filename<<"'."<<std::endl;
      throw std::logic_error("Failed to open a file!");
    }

    #ifdef WITH_COMPRESSION
      boost::iostreams::filtering_ostream out;
      out.push(boost::iostreams::zlib_compressor());
      out.push(fout);
    #else
      auto &out = fout;
    #endif

    out.write(reinterpret_cast<const char*>(&view_height),    sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&view_width),     sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&view_xoff),      sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&view_yoff),      sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&num_data_cells), sizeof(int64_t));
    out.write(reinterpret_cast<const char*>(&no_data),        sizeof(T  ));

    out.write(reinterpret_cast<const char*>(geotransform.data()), 6*sizeof(double));
    std::string::size_type projection_size = projection.size();
    out.write(reinterpret_cast<const char*>(&projection_size), sizeof(std::string::size_type));
    out.write(reinterpret_cast<const char*>(projection.data()), projection.size()*sizeof(const char));

    out.write(reinterpret_cast<const char*>(data.data()), data.size()*sizeof(T));
  }

  ///TODO
  void loadNative(const std::string &filename, bool load_data=true){
    std::ifstream fin(filename, std::ios::in | std::ios::binary);
    assert(fin.good());

    this->filename = filename;
    from_cache    = true;

    #ifdef WITH_COMPRESSION
      boost::iostreams::filtering_istream in;
      in.push(boost::iostreams::zlib_decompressor());
      in.push(fin);
    #else
      auto &in = fin;
    #endif

    in.read(reinterpret_cast<char*>(&view_height),    sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&view_width),     sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&view_xoff),      sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&view_yoff),      sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&num_data_cells), sizeof(int64_t));
    in.read(reinterpret_cast<char*>(&no_data),        sizeof(T  ));
    geotransform.resize(6);
    in.read(reinterpret_cast<char*>(geotransform.data()), 6*sizeof(double));

    std::string::size_type projection_size;
    in.read(reinterpret_cast<char*>(&projection_size), sizeof(std::string::size_type));
    projection.resize(projection_size,' ');
    in.read(reinterpret_cast<char*>(&projection[0]), projection.size()*sizeof(char));

    if(load_data){
      data.resize(view_height*view_width);
      in.read(reinterpret_cast<char*>(data.data()), view_width*view_height*sizeof(T));
    }
  }

 public:
  Array2D(){
    GDALAllRegister();
    view_width   = 0;
    view_height  = 0;
    view_xoff    = 0;
    view_yoff    = 0;
  }

  /**
    @brief Creates a raster of the specified dimensions

    @param[in] width   Width of the raster
    @param[in] height  Height of the raster
    @param[in] val     Initial value of all the raster's cells. Defaults to the
                       Array2D template type's default value
  */
  Array2D(int32_t width, int32_t height, const T& val = T()) : Array2D() {
    resize(width,height,val);
  }

  /**
    @brief Create a raster with the same properties and dimensions as another
           raster. No data is copied between the two.

    @param[in] other   Raster whose properties and dimensions should be copied
    @param[in] val     Initial value of all the raster's cells. Defaults to the
                       Array2D template type's default value
  */
  template<class U>
  Array2D(const Array2D<U> &other, const T& val=T()) : Array2D() {
    view_width   = other.view_width;
    view_height  = other.view_height;
    view_xoff    = other.view_xoff;
    view_yoff    = other.view_yoff;
    geotransform = other.geotransform;
    projection   = other.projection;
    basename     = other.basename;
    resize(other.width(), other.height(), val);
  }

  ///TODO
  Array2D(const std::string &filename, bool native, int32_t xOffset=0, int32_t yOffset=0, int32_t part_width=0, int32_t part_height=0, bool exact=false, bool load_data=true) : Array2D() {
    if(native)
      loadNative(filename, load_data);
    else
      loadGDAL(filename, xOffset, yOffset, part_width, part_height, exact, load_data);
  }

  void setCacheFilename(const std::string &filename){
    this->filename = filename;
  }

  /**
    @brief Caches the raster data and all its properties to disk. Data is then
           purged from RAM.

    @param[in] filename File to save the data to

    @post  Calls to loadData() after this will result in data being loaded from
           the cache.
  */
  void dumpData(){
    saveToCache(filename);
    clear();
  }

  /**
    @brief Loads data from disk into RAM. 

    If dumpData() has been previously called, data is loaded from the cache; 
    otherwise, it is loaded from a GDAL file. No data is loaded if data is
    already present in RAM.
  */
  void loadData() {
    if(!data.empty())
      return; //TODO: Warning?

    if(from_cache){
      loadNative(filename, true);
    } else {
      GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
      if(fin==NULL){
        std::cerr<<"Failed to loadData() into tile from '"<<filename<<"'"<<std::endl;
        throw std::runtime_error("Failed to loadData() into tile.");
      }

      GDALRasterBand *band = fin->GetRasterBand(1);

      data.resize(view_width*view_height);
      auto temp = band->RasterIO( GF_Read, view_xoff, view_yoff, view_width, view_height, data.data(), view_width, view_height, myGDALType(), 0, 0 );
      if(temp!=CE_None)
        throw std::runtime_error("Error reading file with GDAL!");

      GDALClose(fin);
    }
  }

  ///Returns a reference to the internal data array
  T* getData() { return data.data(); }

  ///Number of cells in the DEM
  int64_t size() const { return (int64_t)view_width*(int64_t)view_height; }

  ///Width of the raster
  int32_t width() const { return view_width; }

  ///Height of the raster
  int32_t height() const { return view_height; }

  ///X-Offset of this subregion of whatever raster we loaded from
  int32_t viewXoff() const { return view_xoff; }

  ///Y-Offset of this subregion of whatever raster we loaded from
  int32_t viewYoff() const { return view_yoff; }

  ///Returns TRUE if no data is present in RAM
  bool empty() const { return data.empty(); }

  ///Returns the NoData value of the raster. Cells equal to this value sould
  ///generally not be used in calculations. But note that the isNoData() method
  ///is a much better choice for testing whether a cell is NoData or not.
  T noData() const { return no_data; }

  ///Finds the minimum value of the raster, ignoring NoData cells
  T min() const {
    T minval = std::numeric_limits<T>::max();
    for(auto const x: data){
      if(x==no_data)
        continue;
      minval = std::min(minval,x);
    }
    return minval;
  }

  ///Finds the maximum value of the raster, ignoring NoData cells
  T max() const {
    T maxval = std::numeric_limits<T>::min();
    for(auto const x: data){
      if(x==no_data)
        continue;
      maxval = std::max(maxval,x);
    }
    return maxval;
  }

  /**
    @brief Replace one cell value with another throughout the raster. Can
           operate on NoData cells.

    @param[in] oldval   Value to be replaced
    @param[in] newval   Value to replace 'oldval' with
  */
  void replace(const T oldval, const T newval){
    for(auto &x: data)
      if(x==oldval)
        x = newval;
  }

  /**
    @brief Counts the number of occurrences of a particular value in the raster.
           Can operate on NoData cells.

    @param[in] val   Value to be be counted

    @return The number of times 'val' appears in the raster. Will be 0 if raster
            is not loaded in RAM.
  */
  int64_t countval(const T val) const {
    //TODO: Warn if raster is empty?
    int64_t count=0;
    for(const auto x: data)
      if(x==val)
        count++;
    return count;
  }

  /**
    @brief Convert from index coordinates to x,y coordinates

    @param[in]  i   Index coordinate
    @param[out] x   X-coordinate of i
    @param[out] y   Y-coordinate of i
  */
  void iToxy(const int64_t i, int32_t &x, int32_t &y) const {
    x = i%view_width;
    y = i/view_width;
  }

  /**
    @brief Convert from x,y coordinates to index coordinates

    @param[in]  x   X-coordinate to convert
    @param[in]  y   Y-coordinate to convert

    @return Returns the index coordinate i of (x,y)
  */
  int64_t xyToI(int32_t x, int32_t y) const {
    return (int64_t)y*(int64_t)view_width+(int64_t)x;
  }

  /**
    @brief Given a cell identified by an i-coordinate, return the i-coordinate
           of the neighbour identified by dx,dy

    @param[in]  i   i-coordinate of cell whose neighbour needs to be identified
    @param[in] dx   x-displacement of the neighbour from i
    @param[in] dy   y-displacement of the neighbour from i

    @return i-coordinate of the neighbour. Usually referred to as 'ni'
  */
  int64_t nToI(int64_t i, int32_t dx, int32_t dy) const {
    int32_t x=i%view_width+dx;
    int32_t y=i/view_width+dy;
    if(x<0 || y<0 || x==view_width || y==view_height)
      return -1;
    return xyToI(x,y);
  }

  /**
    @brief Given a cell identified by an i-coordinate, return the i-coordinate
           of the neighbour identified by n

    @param[in]  i   i-coordinate of cell whose neighbour needs to be identified
    @param[in]  n   Neighbour to be identified

    @return i-coordinate of the neighbour. Usually referred to as 'ni'
  */
  int64_t getN(int64_t i, int32_t n) const {
    int32_t x = i%view_width+dx[n];
    int32_t y = i/view_width+dy[n];
    if(x<0 || y<0 || x==view_width || y==view_height)
      return -1;
    return xyToI(x,y);
  }

  /**
    @brief Copies all the properties AND data of another raster into this one

    @param[in]  o   Raster to copy

    @return Returns this raster, with the other raster's data and properties
            copied in.
  */
  template<class U>
  T& operator=(const Array2D<U> &o){
    data = std::vector<T>(o.data.begin(),o.data.end());
    view_height    = o.view_height;
    view_width     = o.view_width;
    view_xoff      = o.view_xoff;
    view_yoff      = o.view_yoff;
    num_data_cells = o.num_data_cells;
    geotransform   = o.geotransform;
    no_data        = (T)o.no_data;
    return *this;
  }

  /**
    @brief Determine if two rasters are equivalent based on dimensions,
           NoData value, and their data
  */
  bool operator==(const Array2D<T> &o){
    if(width()!=o.width() || height()!=o.height())
      return false;
    if(noData()!=o.noData())
      return false;
    return data==o.data;
  }

  /**
    @brief Whether or not a cell is NoData using x,y coordinates

    @param[in]  x   X-coordinate of cell to test
    @param[in]  y   Y-coordinate of cell to test

    @return Returns TRUE if the cell is NoData
  */
  bool isNoData(int32_t x, int32_t y) const {
    assert(0<=x && x<view_width);
    assert(0<=y && y<view_height);
    return data[xyToI(x,y)]==no_data;
  }

  /**
    @brief Whether or not a cell is NoData using i coordinates

    @param[in]  i   i-coordinate of cell to test

    @return Returns TRUE if the cell is NoData
  */
  bool isNoData(int64_t i) const {
    assert(0<=i && i<size());
    return data[i]==no_data;
  }

  /**
    @brief Flips the raster from top to bottom
  */
  void flipVert(){
    for(int32_t y=0;y<view_height/2;y++)
      std::swap_ranges(
        data.begin()+xyToI(0,y),
        data.begin()+xyToI(view_width,y),
        data.begin()+xyToI(0,view_height-1-y)
      );
  }

  /**
    @brief Flips the raster from side-to-side
  */
  void flipHorz(){
    for(int32_t y=0;y<view_height;y++)
      std::reverse(data.begin()+xyToI(0,y),data.begin()+xyToI(view_width,y));
  }

  /**
    @brief Flips the raster about its diagonal axis, like a matrix tranpose.
  */
  void transpose(){
    std::cerr<<"transpose() is an experimental feature."<<std::endl;
    std::vector<T> new_data(view_width*view_height);
    for(int32_t y=0;y<view_height;y++)
    for(int32_t x=0;x<view_width;x++)
      new_data[(int64_t)x*(int64_t)view_height+(int64_t)y] = data[xyToI(x,y)];
    data = new_data;
    std::swap(view_width,view_height);
    //TODO: Offsets?
  }

  /**
    @brief Test whether a cell lies within the boundaries of the raster

    @param[in]  x   X-coordinate of cell to test
    @param[in]  y   Y-coordinate of cell to test

    @return TRUE if cell lies within the raster
  */
  bool inGrid(int32_t x, int32_t y) const {
    return 0<=x && x<view_width && 0<=y && y<view_height;
  }

  /**
    @brief Test whether a cell lies within the boundaries of the raster.

    Obviously this bears some difference from `inGrid(x,y)`.

    @param[in]  i   i-coordinate of cell to test

    @return TRUE if cell lies within the raster
  */
  bool inGrid(int64_t i) const {
    return 0<=i && i<size();
  }

  /**
    @brief Test whether a cell lies on the boundary of the raster

    @param[in]  x   X-coordinate of cell to test
    @param[in]  y   X-coordinate of cell to test

    @return TRUE if cell lies on the raster's boundary
  */
  bool isEdgeCell(int32_t x, int32_t y) const {
    return x==0 || y==0 || x==view_width-1 || y==view_height-1;
  }

  /**
    @brief Test whether a cell lies on the boundary of the raster

    @param[in]  i   i-coordinate of cell to test

    @return TRUE if cell lies on the raster's boundary
  */
  bool isEdgeCell(int64_t i) const {
    int32_t x,y;
    iToxy(i,x,y);
    return isEdgeCell(x,y);
  }

  /**
    @brief Sets the NoData value of the raster

    @param[in]   ndval    Value to change NoData to
  */
  void setNoData(const T &ndval){
    no_data = ndval;
  }

  /**
    @brief Sets all of the raster's cells to 'val'

    @param[in]   val      Value to change the cells to
  */
  void setAll(const T val){
    std::fill(data.begin(),data.end(),val);
  }

  /**
    @brief Resize the raster. Note: this clears all the raster's data.

    @param[in]   width    New width of the raster
    @param[in]   height   New height of the raster
    @param[in]   val      Value to set all the cells to. Defaults to the 
                          raster's template type default value
  */
  void resize(int32_t width, int32_t height, const T& val = T()){
    data.resize(width*height);
    setAll(val);
    view_height = height;
    view_width  = width;
  }

  /*
    @brief Resize a raster to copy another raster's dimensions. Copy properies.

    @param[in]   other    Raster to match sizes with
    @param[in]   val      Value to set all the cells to. Defaults to the
                          raster's template type default value
  */
  template<class U>
  void resize(const Array2D<U> &other, const T& val = T()){
    resize(other.width(), other.height(), val);
    geotransform = other.geotransform;
    projection   = other.projection;
  }

  /**
    @brief Makes a raster larger and retains the raster's old data, similar to resize.

    Note: Using this command requires RAM equal to the sum of the old raster and
    the new raster. The old raster is placed in the upper-left of the new
    raster.

    @param[in] new_width  New width of the raster. Must be >= the old width.
    @param[in] new_height New height of the raster. Must be >= the old height.
    @param[in] val        Value to set the new cells to
  */
  void expand(int32_t new_width, int32_t new_height, const T val){
    if(new_width<view_width)
      throw std::runtime_error("expand(): new_width<view_width");
    if(new_height<view_height)
      throw std::runtime_error("expand(): new_height<view_height");
    
    int32_t old_width  = width();
    int32_t old_height = height();

    std::vector<T> old_data = std::move(data);

    resize(new_width,new_height,val);

    for(int32_t y=0;y<old_height;y++)
    for(int32_t x=0;x<old_width;x++)
      data[y*new_width+x] = old_data[y*old_width+x];
  }

  /**
    @brief Counts the number of cells which are not NoData.
  */
  void countDataCells(){
    num_data_cells = 0;
    for(const auto x: data)
      if(x!=no_data)
        num_data_cells++;
  }

  /**
    @brief Returns the number of cells which are not NoData. May count them.

    @return Returns the number of cells which are not NoData.
  */
  int64_t numDataCells(){
    if(num_data_cells==-1)
      countDataCells();
    return num_data_cells;
  }

  /**
    @brief Returns the number of cells which are not NoData. Does not count them.

    countDataCells() should be call prior to running this method, or the
    non-const version of the method should be used.

    @return Returns the number of cells which are not NoData.
  */
  int64_t numDataCells() const {
    return num_data_cells;
  }

  /**
    @brief Return cell value based on i-coordinate

    @param[in]   i    i-coordinate of cell whose data should be fetched.

    @return The value of the cell identified by 'i'
  */
  T& operator()(int64_t i){
    assert(i>=0);
    assert(i<view_width*view_height);
    return data[i];
  }

  /**
    @brief Return cell value based on i-coordinate

    @param[in]   i    i-coordinate of cell whose data should be fetched.

    @return The value of the cell identified by 'i'
  */
  T operator()(int64_t i) const {
    assert(i>=0);
    assert(i<view_width*view_height);
    return data[i];
  }

  /**
    @brief Return cell value based on x,y coordinates

    @param[in]   x    X-coordinate of cell whose data should be fetched.
    @param[in]   y    Y-coordinate of cell whose data should be fetched.

    @return The value of the cell identified by x,y
  */
  T& operator()(int32_t x, int32_t y){
    assert(x>=0);
    assert(y>=0);
    //std::cerr<<"Width: "<<width()<<" Height: "<<height()<<" x: "<<x<<" y: "<<y<<std::endl;
    assert(x<width());
    assert(y<height());
    return data[xyToI(x,y)];
  }

  /**
    @brief Return cell value based on x,y coordinates

    @param[in]   x    X-coordinate of cell whose data should be fetched.
    @param[in]   y    Y-coordinate of cell whose data should be fetched.

    @return The value of the cell identified by x,y
  */
  T operator()(int32_t x, int32_t y) const {
    assert(x>=0);
    assert(y>=0);
    assert(x<width());
    assert(y<height());
    return data[xyToI(x,y)];
  }

  /**
    @brief Returns a copy of the top row of the raster

    @return A vector containing a copy of the top row of the raster
  */
  std::vector<T> topRow() const {  
    return getRowData(0);  
  }

  /**
    @brief Returns a copy of the bottom row of the raster

    @return A vector containing a copy of the bottom row of the raster
  */
  std::vector<T> bottomRow() const {
    return getRowData(view_height-1);
  }

  /**
    @brief Returns a copy of the left column of the raster

    Top to bottom is reoriented as left to right.

    @return A vector containing a copy of the left column of the raster
  */
  std::vector<T> leftColumn() const { 
    return getColData(0); 
  }

  /**
    @brief Returns a copy of the right column of the raster

    Top to bottom is reoriented as left to right.

    @return A vector containing a copy of the right column of the raster
  */
  std::vector<T> rightColumn() const { 
    return getColData(view_width-1);
  }

  /**
    @brief Sets an entire row of a raster to a given value.

    @param[in]   y    The row to be set
    @param[in] val    The value to set the row to
  */
  void setRow(int32_t y, const T &val){
    std::fill(data.begin()+xyToI(y,0),data.begin()+xyToI(y,view_width),val);
  }

  /**
    @brief Sets an entire column of a raster to a given value.

    @param[in]   x    The column to be set
    @param[in] val    The value to set the column to
  */
  void setCol(int32_t x, const T &val){
    for(int32_t y=0;y<view_height;y++)
      data[xyToI(x,y)] = val;
  }

  /**
    @brief Returns a copy of an arbitrary row of the raster

    @param[in]   y    The row to retrieve

    @return A vector containing a copy of the selected row
  */
  std::vector<T> getRowData(int32_t y) const {
    return std::vector<T>(data.begin()+xyToI(0,y),data.begin()+xyToI(view_width,y));
  }

  /**
    @brief Returns a copy of an arbitrary column of the raster

    @param[in]   x    The column to retrieve

    @return A vector containing a copy of the selected column
  */
  std::vector<T> getColData(int32_t x) const {
    std::vector<T> temp(view_height);
    for(int32_t y=0;y<view_height;y++)
      temp[y]=data[xyToI(x,y)];
    return temp;
  }

  ///Clears all raster data from RAM
  void clear(){
    data.clear();
    data.shrink_to_fit();
  }

  /**
    @brief Copies the geotransform, projection, and basename of another raster

    @param[in]    other    Raster to copy from
  */
  template<class U>
  void templateCopy(const Array2D<U> &other){
    geotransform = other.geotransform;
    projection   = other.projection;
    basename     = other.basename;
  }

  void saveGDAL(const std::string &filename, int32_t xoffset, int32_t yoffset){
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if(poDriver==NULL){
      std::cerr<<"Could not open GDAL driver!"<<std::endl;
      throw std::runtime_error("Could not open GDAL driver!");
    }
    GDALDataset *fout    = poDriver->Create(filename.c_str(), width(), height(), 1, myGDALType(), NULL);
    if(fout==NULL){
      std::cerr<<"Could not open file '"<<filename<<"' for GDAL save!"<<std::endl;
      throw std::runtime_error("Could not open file for GDAL save!");
    }

    GDALRasterBand *oband = fout->GetRasterBand(1);
    oband->SetNoDataValue(no_data);

    //The geotransform maps each grid cell to a point in an affine-transformed
    //projection of the actual terrain. The geostransform is specified as follows:
    //    Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)
    //    Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)
    //In case of north up images, the GT(2) and GT(4) coefficients are zero, and
    //the GT(1) is pixel width, and GT(5) is pixel height. The (GT(0),GT(3))
    //position is the top left corner of the top left pixel of the raster.

    auto out_geotransform = geotransform;

    if(out_geotransform.size()!=6){
      std::cerr<<"Geotransform of output is not the right size. Found "<<out_geotransform.size()<<" expected 6."<<std::endl;
      throw std::runtime_error("saveGDAL(): Invalid output geotransform.");
    }

    //We shift the top-left pixel of hte image eastward to the appropriate
    //coordinate
    out_geotransform[0] += xoffset*geotransform[1];

    //We shift the top-left pixel of the image southward to the appropriate
    //coordinate
    out_geotransform[3] += yoffset*geotransform[5];

    fout->SetGeoTransform(out_geotransform.data());
    fout->SetProjection(projection.c_str());

    #ifdef DEBUG
      std::cerr<<"Filename: "<<std::setw(20)<<filename<<" Xoffset: "<<std::setw(6)<<xoffset<<" Yoffset: "<<std::setw(6)<<yoffset<<" Geotrans0: "<<std::setw(10)<<std::setprecision(10)<<std::fixed<<geotransform[0]<<" Geotrans3: "<<std::setw(10)<<std::setprecision(10)<<std::fixed<<geotransform[3]<< std::endl;
    #endif

    auto temp = oband->RasterIO(GF_Write, 0, 0, view_width, view_height, data.data(), view_width, view_height, myGDALType(), 0, 0);
    if(temp!=CE_None)
      std::cerr<<"Error writing file! Continuing in the hopes that some work can be salvaged."<<std::endl;

    GDALClose(fout);
  }

  /**
    @brief Output a square of cells useful for determining raster orientation.

    This method prints out a square block of cells whose upper-left corner is
    the (integer-division) center of the raster.

    Stamps are only shown if the SHOW_STAMPS preprocessor variable is set.

    Since algorithms may have to flip rasters horizontally or vertically before
    manipulating them, it is important that all algorithms work on data in the
    same orientation. This method, used in testing, helps a user ensure that 
    their algorithm is orientating data correctly.

    @param[in]  size   Output stamp will be size x size
  */
  void printStamp(int32_t size, std::string msg="") const {
    #ifdef SHOW_STAMPS
      const int32_t sx = width()/2;
      const int32_t sy = height()/2;

      if(msg.size()>0)
        std::cerr<<msg<<std::endl;
      std::cerr<<"Stamp for basename='"<<basename
               <<"', filename='"<<filename
               <<"', dtype="<<GDALGetDataTypeName(myGDALType())
               <<" at "<<sx<<","<<sy<<"\n";

      const int32_t sxmax = std::min(width(), sx+size);
      const int32_t symax = std::min(height(),sy+size);

      for(int32_t y=sy;y<symax;y++){
        for(int32_t x=sx;x<sxmax;x++)
          std::cerr<<std::setw(5)<<std::setprecision(3)<<(int)data[xyToI(x,y)]<<" ";
        std::cerr<<"\n";
      }
    #endif
  }


  /**
    @brief Prints a square of cells centered at x,y. Useful for debugging.

    @param[in]  radius   Output stamp will be 2*radius x 2*radius
    @param[in]       x   X-coordinate of block center
    @param[in]       y   Y-coordinate of block center
    @parma[in]     msg   Optional message to print above the block
  */
  void printBlock(const int32_t radius, const int32_t x0, const int32_t y0, bool color=false, const std::string msg="") const {
    if(msg)
      std::cerr<<msg<<std::endl;

    int32_t xmin = std::max(0,x0-size);
    int32_t ymin = std::max(0,y0-size);
    int32_t xmax = std::min(width(),x0+size);
    int32_t ymax = std::min(height(),y0+size);

    for(int32_t y=ymin;y<ymax;y++){
      for(int32_t x=xmin;x<xmax;x++){
        if(color && x==x0 && y==y0)
          std::cerr<<"\033[92m";
        std::cerr<<std::setw(5)<<(int)data[xyToI(x,y)]<<" ";
        if(color && x==x0 && y==y0)
          std::cerr<<"\033[39m";
      }
      std::cerr<<std::endl;
    }
  }

  /**
    @brief Get the area of an individual cell in square projection units

    @return The area of the cell in square projection units
  */
  double getCellArea() const {
    return geotransform[1]*geotransform[5];
  }
};

#endif