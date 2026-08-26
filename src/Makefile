CXX = g++
CXXFLAGS = -fopenmp -std=c++20 -O3
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
MAKEFLAGS += -j8

ClonalFrame2: $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) ClonalFrame2
