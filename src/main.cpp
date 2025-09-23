#include <iostream>

#include "Move.h"
#include "Param.h"
#include "Alignment.h"

int main() {
    auto *p = new Param();
    auto *m = new Move(0.5);
    auto *a = new Alignment("test.fasta");
    p->set_alpha(1);
    p->set_beta(2);
    p->set_gamma(3);
    m->action(p);
    p->print();
    delete p;
    delete m;
    delete a;
    return 0;
}