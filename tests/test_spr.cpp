#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

//This testing main function reads the alignment in tmp.fasta and the tree in tmp.nwk, performs parsimonious reconstruction, outputs it, and returns the likelihood
int main(int argc, const char *argv[]) {
    int n=argc>1?std::stoi(argv[1]):4;
    int niter=argc>2?std::stoi(argv[2]):1000;
    auto t = Tree::random_tree(n,0.01);
    for (int i=0;i<niter;i++) {
        std::println("{}",t->to_string());
        for (int j=0;j<10;j++) t->rSPR();
    }
}