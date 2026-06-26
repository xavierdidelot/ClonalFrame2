#include "../src/Param.h"
#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

int main() {
    setSeed(0);
    //omp_set_num_threads(1);
    verbose=true;
    profiling=true;
    auto t = Tree::random_tree(10,0.01,0.01);
    t->simulate_sequences(10000);
    auto p = new Param(t);
    p->mcmc(10,10,1);
    t=p->get_tree();
    delete(p);
    delete(t);
}