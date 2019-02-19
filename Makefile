PROG  := parallelquicksort
PROG2 := makequicksortdata

SRCS  := parallelquicksort.cpp
SRCS2 := makequicksortdata.cpp

OBJS  = parallelquicksort.o
OBJS2 = makequicksortdata.o

DEPS  = parallelquicksort.d
DEPS2 = makequicksortdata.d

VPATH  = src/parallelquicksort src/makequicksortdata
CXX = icpc
CXXFLAGS = -Wall -Wextra -O3 -xHOST -ipo -pipe -std=c++17 -fopenmp -I/home/dc1394/oss/parallelstl-20181204/include
LDFLAGS = -ltbb \
		  -L/home/dc1394/oss/boost_1_69_0/stage/icc/lib -lboost_filesystem -lboost_system -lboost_thread

all: $(PROG) $(PROG2) ;
#rm -f $(OBJS) $(DEPS)

-include $(DEPS)

$(PROG): $(OBJS)
		$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^

$(PROG2): $(OBJS2)
		$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
		$(CXX) $(CXXFLAGS) -c -MMD -MP -DDEBUG $<

clean:
		rm -f $(PROG) $(OBJS) $(DEPS)
		rm -f $(PROG2) $(OBJS2) $(DEPS2)
