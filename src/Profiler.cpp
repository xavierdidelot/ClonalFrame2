//
// Created by Didelot, Xavier on 05/05/2026.
//

#include "Profiler.h"

Profiler Profiler::prof;
thread_local std::unordered_map<const char*, Profiler::Clock::time_point> Profiler::tl_start_times;

static Profiler::Entry& getOrInsert(
    std::unordered_map<const char*, Profiler::Entry>& entries,
    std::mutex& mtx,
    const char* name)
{
    {
        auto it = entries.find(name);
        if (it != entries.end()) return it->second;
    }
    std::lock_guard<std::mutex> lk(mtx);
    auto [it, _] = entries.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(name),
        std::forward_as_tuple());
    it->second.name = name;
    return it->second;
}

void Profiler::start(const char* name) {
    tl_start_times[name] = Clock::now();
    getOrInsert(entries, entries_mutex, name);
}

void Profiler::stop(const char* name) {
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(Clock::now() - tl_start_times[name]).count();
    Entry& e = getOrInsert(entries, entries_mutex, name);
    e.duration.fetch_add(us, std::memory_order_relaxed);
    e.calls.fetch_add(1,   std::memory_order_relaxed);
}

Profiler::Profiler() {
    total_start = Clock::now();
}

Profiler::~Profiler() {
    if (profiling) prof.print();
}

void Profiler::print() {
    std::vector<Entry*> sorted;
    sorted.reserve(entries.size());
    for (auto& [_, e] : entries) sorted.push_back(&e);

    std::sort(sorted.begin(), sorted.end(),
        [](const Entry* a, const Entry* b) {
            return a->duration.load() < b->duration.load();
        });

    long long total = std::chrono::duration_cast<std::chrono::microseconds>(
        Clock::now() - total_start).count();

    std::cout<<"\n=== Profiler Results ===\n";
    std::cout << std::left
          << std::setw(40) << "Function" << ' '
          << std::setw(12) << "Calls" << ' '
          << std::setw(16) << "Total (ms)" << ' '
          << std::setw(16) << "Percentage" << '\n';
    std::cout << std::string(84, '-') << '\n';

    for (const Entry* e : sorted) {
        std::cout << std::left << std::setw(40) << e->name << ' '
                  << std::setw(12) << e->calls.load() << ' '
                  << std::fixed << std::setprecision(2)
                  << std::setw(16) << e->duration.load() / 1000.0 << ' '
                  << std::setw(16) << static_cast<double>(e->duration.load()) / total * 100.0
                  << '\n';
    }
}