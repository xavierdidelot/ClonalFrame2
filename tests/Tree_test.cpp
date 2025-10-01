#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>

#include <string>
#include <fstream>

using namespace std;

int main() {
    auto *t = new Tree(10);
    assert(t->get_n()==10);
    t->write("tmp.nwk");
    //system("open tmp.nwk");
    remove("tmp.nwk");
    t->add_mutations(100,10);
    assert(t->check_structure());
    auto * a=new Alignment(t);
    a->write("tmp.fasta");
    delete(a);
    delete(t);
}