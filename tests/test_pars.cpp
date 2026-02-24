#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

using namespace std;

//This testing main function reads the alignment in tmp.fasta and the tree in tmp.nwk, performs parsimonious reconstruction and outputs the number of substitutions
int main() {
    auto a = new Alignment("tmp.fasta");
    auto t = new Tree("tmp.nwk");
    t->set_seqs(a);
    int s1=t->ancestral_reconstruction();
    int s2=t->count_subs();
    cout<<s1<<","<<s2<<endl;
    delete(t);
    delete(a);
}