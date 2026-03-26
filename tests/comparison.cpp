#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

//This script compares the likelihood estimated using several approximations
//Arguments are n, l, mu, niter, eps
int main(int argc, const char *argv[]) {
    int n=argc>1?std::stoi(argv[1]):20;
    int l=argc>2?std::stoi(argv[2]):10000;
    double mu=argc>3?std::stof(argv[3]):0.01;
    int niter=argc>4?std::stoi(argv[4]):100;
    double eps=argc>5?std::stof(argv[5]):0.01;
    auto t = Tree::random_tree(n,mu);
    t->simulate_sequences(l);
    auto a=Alignment::make_from_tree(t,false,false);
    std::ofstream file;
    //file.open("tmp.nwk");
    std::println(file,"{}",t->to_string());
    //a->write("tmp.fasta",true);
    for (int i=0;i<niter;i++) {
        auto t2=Tree::noise_tree(t,eps);
        std::println(file,"{}",t2->to_string());
        t2->set_sequences(a);
        double r1=t2->likelihood();
        double r2=t2->jointml();
        t2->parsimony();
        double r3=t2->likelihood_anc();
        std::println("{:.10e},{:.10e},{:.10e}", r1, r2, r3);
    }
    file.close();
}