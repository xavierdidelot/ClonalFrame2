#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

using namespace std;

//This testing main function reads the alignment in tmp.fasta and the tree in tmp.nwk, performs parsimonious reconstruction, outputs it, and returns the likelihood
int main() {
    auto a = Alignment::read_fasta("tmp.fasta");
    auto t = Tree::read_newick("tmp.nwk");
    t->set_sequences(a);
    double l = t->likelihood();
    t->parsimony();
    a=Alignment::make_from_tree(t,true,false);
    t->write("tmp2.nwk");
    a->write("tmp2.fasta",false);
    double s3=t->likelihood_anc();
    cout<<std::fixed << std::setprecision(6)<<l<<","<<std::fixed << std::setprecision(6)<<s3<<endl;
}