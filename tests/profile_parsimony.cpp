#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    auto t = Tree::random_tree(1000,1e-4,1e-4);
    t->simulate_sequences(10000);
    int s1=0,s2=0;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    myprintln("Simulation time: {} seconds",elapsed.count());

    start = std::chrono::high_resolution_clock::now();
    for (int i=0;i<100;i++) {
        s1=t->parsimony();
    }
    elapsed = std::chrono::high_resolution_clock::now() - start;
    myprintln("Parsimony time: {} seconds",elapsed.count());

    start = std::chrono::high_resolution_clock::now();
    for (int i=0;i<100;i++) {
        s2=t->parsimony_recursive();
    }
    elapsed = std::chrono::high_resolution_clock::now() - start;
    myprintln("Parsimony recursive time: {} seconds",elapsed.count());

    myprintln("{} {}",s1,s2);
    delete(t);
}