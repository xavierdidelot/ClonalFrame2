#ifndef CLONALFRAME2_GLOBAL_H
#define CLONALFRAME2_GLOBAL_H
#ifdef _OPENMP
#include <omp.h>
#endif
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "Profiler.h"

inline bool verbose=false; //Verbose mode

inline std::ofstream out; //Output file

inline std::string vtime()
{
    auto time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::ostringstream out; out << std::put_time(std::localtime(&time), "[%H:%M:%S] ");
    return out.str();
}

inline std::mt19937 rng(std::random_device{}()); // Random number generator, using a non-deterministic seed by default

//Function to set the seed of the random number generator
inline void setSeed(unsigned int seed) {
    rng.seed(seed);
}

inline int max_available_threads() {
    #ifdef _OPENMP
        return omp_get_max_threads();
    #else
        return 1;
    #endif
}

inline int n_cores = max_available_threads();//Number of cores to use, defaulting to OpenMP default or OpenMP is available or zero otherwise

inline void set_num_threads(int count) {
    #ifdef _OPENMP
        omp_set_num_threads(count);
        n_cores=count;
    #else
       std::cout<<"Ignoring requested number of cores since OpenMP is not available"<<std::endl;
    #endif
}


inline bool profiling=false;//Profiling mode

inline const size_t maxskip=1000;

///dgamma functions, returning log(pdf)
inline double dgamma_shape_scale(double x,double shape,double scale) {
    if (x==0) return 0.0;//this is not generally correct but will cancel out in calculations
    return (shape - 1.0) * log(x) - x / scale - shape * log(scale) - lgamma(shape);
}

inline double dgamma_mean_sd(double x,double mean,double sd) {
    double shape=mean*mean/sd/sd;
    double scale=sd*sd/mean;
    if (x==0) return 0.0;//this is not generally correct but will cancel out in calculations
    return (shape - 1.0) * log(x) - x / scale - shape * log(scale) - lgamma(shape);
}

#endif
