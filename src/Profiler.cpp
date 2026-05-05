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
    entries[name].total_us += us;
    entries[name].calls    += 1;
}

Profiler::~Profiler() {
    this->print();
}

void Profiler::print() {
    if (!profiling) return;
    std::vector<std::pair<const char*, Entry>> sorted(entries.begin(), entries.end());
    std::sort(sorted.begin(), sorted.end(),
        [](const auto& a, const auto& b) {
            return a.second.total_us > b.second.total_us;
        });

    std::println("\n=== Profiler Results ===");
    std::println("{:<40} {:<12} {:<16} {:<16}", "Function", "Calls", "Total (ms)", "Avg (µs)");
    std::println("{}", std::string(84, '-'));
    for (const auto& [name, e] : sorted) {
        std::println("{:<40} {:<12} {:<16.2f} {:<16.2f}",
            name,
            e.calls,
            e.total_us / 1000.0,
            static_cast<double>(e.total_us) / e.calls);
    }
}