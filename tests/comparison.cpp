#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

//This script compares the likelihood estimated using several approximations
//Arguments are n, l, mu, niter, seed
int main(int argc, const char *argv[]) {
    if (argc<6) return 0;
    int n=std::stoi(argv[1]);
    int l=std::stoi(argv[2]);
    double avg_bra_len=std::stof(argv[3]);
    int niter=std::stoi(argv[4]);
    int seed=std::stoi(argv[5]);
    setSeed(seed);
    auto t = Tree::random_tree(n,avg_bra_len,avg_bra_len);
    t->simulate_sequences(l);
    auto a=Alignment::make_from_tree(t,false);
    //a->write("tmp.fasta",true);
    for (int i=0;i<niter;i++) {
        Tree * t2;
        if (i==0) t2=Tree::copy_tree(t); else {
            t2=Tree::noise_tree(t,avg_bra_len);
            t2->set_sequences(a);
        }
        std::println("{}",t2->to_string());
        double r1=t2->likelihood();
        double r2=t2->jointml();
        t2->parsimony();
        double r3=t2->likelihood_anc();
        std::println("{:.10e}\t{:.10e}\t{:.10e}", r1, r2, r3);
        delete t2;
    }
}