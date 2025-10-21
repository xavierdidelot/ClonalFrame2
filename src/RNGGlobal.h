#ifndef RNG_GLOBAL_H
#define RNG_GLOBAL_H
#include <random>

using namespace std;

inline bool verbose=false; // Defined as global variable

inline mt19937 rng(random_device{}()); // Global RNG, non-deterministic seed by default

inline void setSeed(unsigned int seed) {
    rng.seed(seed);
}

#endif