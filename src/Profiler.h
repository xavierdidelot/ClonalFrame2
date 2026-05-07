//
// Created by Didelot, Xavier on 05/05/2026.
//

#ifndef CLONALFRAME2_PROFILER_H
#define CLONALFRAME2_PROFILER_H


#include "RNGGlobal.h"
#include <chrono>
#include <unordered_map>
#include <atomic>
#include <mutex>

class Profiler {
public:
    using Clock = std::chrono::high_resolution_clock;
    struct Entry {
        std::atomic<long long> duration{0};
        std::atomic<long long> calls{0};
        const char* name = nullptr;
    };
    static Profiler prof;
    void start(const char* name);
    void stop(const char* name);
    void print();

private:
    std::unordered_map<const char*, Entry> entries;
    std::mutex entries_mutex;
    Clock::time_point total_start;  // plain member, safe at static init time
    static thread_local std::unordered_map<const char*, Clock::time_point> tl_start_times;
    Profiler();
    ~Profiler();
};

#endif //CLONALFRAME2_PROFILER_H