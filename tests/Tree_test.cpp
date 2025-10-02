#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>

#include <string>
#include <fstream>

using namespace std;

int main() {
    auto * t = new Tree(10);
    assert(t->get_n()==10);
    assert(t->check_structure());
    //t = new Tree(10,1);
    assert(t->get_n()==10);
    assert(t->check_structure());
    t->write("tmp.nwk");
    //system("open tmp.nwk");
    //remove("tmp.nwk");
    t->add_mutations(10000,200);
    auto * a=new Alignment(t);
    a->write("tmp.fasta");
    auto * t2=new Tree(a);
    assert(t2->get_n()==10);
    assert(t2->check_structure());
    t2->write("tmp2.nwk");
    delete(a);
    delete(t);
    delete(t2);
}