#ifndef RNG_GLOBAL_H
#define RNG_GLOBAL_H
#include <random>
#include <print>
#include <chrono>

inline bool verbose=false; // Defined as global variable

template<typename... Args>
void vprint(std::format_string<Args...> fmt, Args&&... args)
{
    auto now = std::chrono::system_clock::now();
    std::println("[{:%H:%M:%S}]  {}", now, std::format(fmt, std::forward<Args>(args)...));
}

inline std::mt19937 rng(std::random_device{}()); // Global RNG, non-deterministic seed by default

inline void setSeed(unsigned int seed) {
    rng.seed(seed);
}

#endif