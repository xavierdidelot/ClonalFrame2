#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>

#include <string>
#include <fstream>

using namespace std;

int main() {
    auto *t = new Tree(10);
    assert(t->get_n()==10);
    assert(t->check_structure());
    //assert(t->get_sum_bra_len()==8);
    t = new Tree(10,1);
    assert(t->get_n()==10);
    assert(t->check_structure());
    t->write("tmp.nwk");
    //system("open tmp.nwk");
    //remove("tmp.nwk");
    t->add_mutations(1000,200);
    auto * a=new Alignment(t);
    a->write("tmp.fasta");
    delete(a);
    delete(t);
}