#ifndef RNG_GLOBAL_H
#define RNG_GLOBAL_H
#include <omp.h>
#include <random>
#include <print>
#include <chrono>
#include "Profiler.h"

inline bool verbose=false; //Verbose mode, defined as global variable

template<typename... Args>
void vprint(std::format_string<Args...> fmt, Args&&... args)
{
    auto now = std::chrono::system_clock::now();
    auto ms = floor<std::chrono::milliseconds>(now);
    std::println("[{:%H:%M:%S}]  {}", ms, std::format(fmt, std::forward<Args>(args)...));
}

inline std::mt19937 rng(std::random_device{}()); // Global RNG, non-deterministic seed by default

inline void setSeed(unsigned int seed) {
    rng.seed(seed);
}

inline int n_cores = omp_get_max_threads();//Number of cores to use, defined as a global variable and defaulting to OpenMP default

inline bool profiling=false;//Profiling move, defined as global variable

inline const size_t maxskip=1000;

#endif