#include "../src/Param.h"
#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

using namespace std;

int main() {
    auto t = new Tree(10);
    t->randomise(1.0);
    t->add_mutations(10000,10);
    auto p = new Param(t);
    p->mcmc(10,10,1);
    delete(p);
    delete(t);
}