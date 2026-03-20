#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

//This testing main function reads the alignment in tmp.fasta, creates a NJ tree and writes it into nj.nwk
int main() {
    auto a = Alignment::read_fasta("tmp.fasta");
    auto mat=a->distmat();
    auto t = Tree::nj(mat,a->get_names());
    t->write("nj.nwk");
    delete(mat);
    delete(t);
    delete(a);
}