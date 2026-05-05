//
// Created by Didelot, Xavier on 05/05/2026.
//

#ifndef CLONALFRAME2_PROFILER_H
#define CLONALFRAME2_PROFILER_H


#include "RNGGlobal.h"
#include <chrono>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <print>
#include <vector>


class Profiler {
public:
    static Profiler prof;
    Profiler() = default;
    ~Profiler();
    void print();
    void start(const char* name);
    void stop(const char* name);

private:
    using Clock = std::chrono::high_resolution_clock;
    struct Entry {
        long long total_us = 0;
        long long calls    = 0;
        Clock::time_point last_started;
    };
    std::unordered_map<const char*, Entry> entries;
};

#endif //CLONALFRAME2_PROFILER_H