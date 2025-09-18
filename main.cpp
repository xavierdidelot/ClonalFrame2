#include <iostream>

#include "Move.h"
#include "Param.h"

int main() {
    auto lang = "C++";
    auto *p = new Param();
    auto *m = new Move(0.5);
    p->set_alpha(1);
    p->set_beta(2);
    p->set_gamma(3);
    std::cout << "Hello and welcome to " << lang << "!\n";

    for (int i = 1; i <= 5; i++) {
        std::cout << "i = " << i << std::endl;
    }

    m->action(p);
    p->print();
    delete p;
    delete m;
    return 0;
}