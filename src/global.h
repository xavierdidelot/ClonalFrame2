#ifndef GLOBAL_H
#define GLOBAL_H
#include <omp.h>
#include <random>
#include <print>
#include <chrono>
#include <fstream>
#include <sstream>
#include <iostream>
#include "Profiler.h"

inline bool verbose=false; //Verbose mode

inline std::ofstream out; //Output file

//Function to write a message
template<typename... Args>
void vprint(std::format_string<Args...> fmt, Args&&... args)
{
    auto now = std::chrono::system_clock::now();
    auto ms = floor<std::chrono::milliseconds>(now);
    std::println("[{:%H:%M:%S}]  {}", ms, std::format(fmt, std::forward<Args>(args)...));
}

inline std::mt19937 rng(std::random_device{}()); // Random number generator, using a non-deterministic seed by default

//Function to set the seed of the random number generator
inline void setSeed(unsigned int seed) {
    rng.seed(seed);
}

inline int n_cores = omp_get_max_threads();//Number of cores to use, defaulting to OpenMP default

inline bool profiling=false;//Profiling mode

inline const size_t maxskip=1000;

///dgamma functions, returning log(pdf)
inline double dgamma_shape_scale(double x,double shape,double scale)
{return (shape - 1.0) * log(x) - x / scale - shape * log(scale) - lgamma(shape);}

inline double dgamma_mean_sd(double x,double mean,double sd) {
    double shape=mean*mean/sd/sd;
    double scale=sd*sd/mean;
    return (shape - 1.0) * log(x) - x / scale - shape * log(scale) - lgamma(shape);
}

#endif