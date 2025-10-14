#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>

#include <string>

using namespace std;

int main() {
    auto * t = new Tree(10);
    int n=t->get_n();assert(n==10);
    bool b=t->check_structure();assert(b);
    t->randomise(1.0);
    b=t->check_structure();assert(b);
    t->write("tmp.nwk");
    t->add_mutations(1000,200);
    auto * a=new Alignment(t);
    a->write("tmp.fasta");
    t->nj(a);
    b=t->check_structure();assert(b);
    t->write("tmp2.nwk");
    t->set_seqs(a);
    t->ancestral_reconstruction();
    delete(t);
    delete(a);
}