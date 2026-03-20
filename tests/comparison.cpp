#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

//This script compares the likelihood estimated using several approximations
int main() {
    auto t = Tree::random_tree(20,0.01);
    t->simulate_sequences(10000);
    auto a=Alignment::make_from_tree(t,false,false);//is the last arg always false? if so remove it
    for (int i=0;i<100;i++) {
        auto t2=Tree::noise_tree(t,0.01);
        t2->set_sequences(a);
        double r1=t2->likelihood();
        t2->parsimony();
        double r2=t2->likelihood_anc();
        double r3=t2->jointml();
        std::println("{:.6e},{:.6e},{:.6e}", r1, r2, r3);
    }
}