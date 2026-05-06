//
// Created by Didelot, Xavier on 05/05/2026.
//

#include "Profiler.h"
#include "RNGGlobal.h"

Profiler Profiler::prof;

void Profiler::start(const char* name) {
    entries[name].last_started=Clock::now();
}

void Profiler::stop(const char* name) {
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(
        Clock::now() - entries[name].last_started).count();
    entries[name].duration += us;
    entries[name].calls    += 1;
    entries[name].name = name;
}

Profiler::Profiler() {
    start("Total");
}

Profiler::~Profiler() {
    stop("Total");
    this->print();
}

void Profiler::print() {
    if (!profiling) return;
    std::vector<std::pair<const char*, Entry>> sorted(entries.begin(), entries.end());
    std::sort(sorted.begin(), sorted.end(),
        [](const auto& a, const auto& b) {
            return a.second.duration < b.second.duration;
        });

    std::println("\n=== Profiler Results ===");
    std::println("{:<40} {:<12} {:<16} {:<16}", "Function", "Calls", "Total (ms)", "Percentage");
    std::println("{}", std::string(84, '-'));
    for (const auto& [duration, e] : sorted) {
        std::println("{:<40} {:<12} {:<16.2f} {:<16.2f}",
            e.name,
            e.calls,
            e.duration / 1000.0,
            static_cast<double>(e.duration) / entries["Total"].duration*100.0);
    }
}