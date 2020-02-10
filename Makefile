PROG  := parallelquicksort
PROG2 := makequicksortdata

SRCS  := parallelquicksort.cpp
SRCS2 := makequicksortdata.cpp

OBJS  = parallelquicksort.o
OBJS2 = makequicksortdata.o

DEPS  = parallelquicksort.d
DEPS2 = makequicksortdata.d

VPATH  = src/parallelquicksort src/makequicksortdata
CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -mtune=native -march=native -pipe -std=c++17 -fopenmp
LDFLAGS = -ltbb \
		  -L/home/dc1394/oss/boost_1_72_0/stage/gcc/lib -lboost_filesystem -lboost_system -lboost_thread

all: $(PROG) $(PROG2) ;
#rm -f $(OBJS) $(DEPS)

-include $(DEPS)

$(PROG): $(OBJS)
		$(CXX) $^ $(CXXFLAGS) $(LDFLAGS) -o $@

$(PROG2): $(OBJS2)
		$(CXX) $^ $(CXXFLAGS) $(LDFLAGS) -o $@

%.o: %.cpp
		$(CXX) $(CXXFLAGS) -c -MMD -MP -DDEBUG $<

clean:
		rm -f $(PROG) $(OBJS) $(DEPS)
		rm -f $(PROG2) $(OBJS2) $(DEPS2)
