#include "../src/Tree.h"
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
    assert(t->check_structure());
    delete(t);
}