#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

//This testing main function reads the alignment in tmp.fasta and the tree in tmp.nwk, performs parsimonious reconstruction, outputs it, and returns the likelihood
int main() {
    auto a = Alignment::read_fasta("tmp.fasta");
    auto t = Tree::read_newick("tmp.nwk");
    t->set_sequences(a);
    double l = t->likelihood();
    t->parsimony();
    a=Alignment::make_from_tree(t,true);
    t->write("tmp2.nwk");
    a->write("tmp2.fasta",true);
    a->write("tmp3.fasta",false);
    double s3=t->likelihood_anc();
    double s4=t->likelihood_anc_rec(0.0,202,0.055);
    double s5=t->likelihood_anc_rec(1.1,202,0.055);
    std::cout<<std::setprecision(17)<<l<<","<<s3<<","<<s4<<","<<s5<<std::endl;
}