#ifndef RNG_GLOBAL_H
#define RNG_GLOBAL_H
#include <random>

inline bool verbose=false; // Defined as global variable

inline std::mt19937 rng(std::random_device{}()); // Global RNG, non-deterministic seed by default

inline void setSeed(unsigned int seed) {
    rng.seed(seed);
}

#endif