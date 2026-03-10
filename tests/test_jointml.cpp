#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

using namespace std;

//This testing main function reads the alignment in tmp.fasta and the tree in tmp.nwk, performs parsimonious reconstruction and outputs the number of substitutions
int main() {
    auto a = Alignment::read_fasta("tmp.fasta");
    auto t = Tree::read_newick("tmp.nwk");
    t->set_sequences(a);
    double s1=t->jointml();
    double s2=t->likelihood_anc();
    cout<<std::fixed<<std::setprecision(6)<<s1<<","<<std::fixed<<std::setprecision(6)<<s2<<endl;
    delete(t);
    delete(a);
}