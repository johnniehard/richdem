export CXX=g++
export GDAL_LIBS=`gdal-config --libs`
export GDAL_CFLAGS=`gdal-config --cflags` -DUSEGDAL
RICHDEM_GIT_HASH=`git rev-parse HEAD`
RICHDEM_COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`
export CXXFLAGS=$(GDAL_CFLAGS) --std=c++11 -I../../include -I. -Wall -Wno-unknown-pragmas -DRICHDEM_GIT_HASH="\"$(RICHDEM_GIT_HASH)\"" -DRICHDEM_COMPILE_TIME="\"$(RICHDEM_COMPILE_TIME)\""
export OPT_FLAGS=-O3 -g
export DEBUG_FLAGS=-g -O0 -fno-omit-frame-pointer -DGLIBCXX_DEBUG

#-DNDEBUG -DSHOW_STAMPS

compile: main.cpp
	$(CXX) $(OPT_FLAGS) $(CXXFLAGS) -o breaching main.cpp $(GDAL_LIBS)

debug: main.cpp
	$(CXX) $(DEBUG_FLAGS) $(CXXFLAGS) -o breaching main.cpp $(GDAL_LIBS)