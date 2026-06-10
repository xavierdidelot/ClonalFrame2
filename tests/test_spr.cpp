#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

//This testing main function performs the rSPR move many times
int main(int argc, const char *argv[]) {
    int n=argc>1?std::stoi(argv[1]):4;
    int niter=argc>2?std::stoi(argv[2]):1000;
    int seed=argc>3?std::stoi(argv[3]):1;
    setSeed(seed);
    auto t = Tree::random_tree(n,0.01);
    t->simulate_sequences(0);
    for (int i=0;i<niter;i++) {
        std::println("{}",t->to_string());
        string s;
        for (int j=0;j<10;j++) t->rSPR(&s);
    }
}