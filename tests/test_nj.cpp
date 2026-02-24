#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

using namespace std;

//This testing main function reads the alignment in tmp.fasta, creates a NJ tree and writes it into nj.nwk
int main() {
    auto a = Alignment::read_fasta("tmp.fasta");
    auto t = Tree::nj(a);
    t->write("nj.nwk");
    delete(t);
    delete(a);
}