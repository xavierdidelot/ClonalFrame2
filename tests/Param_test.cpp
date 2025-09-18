#include "../src/Param.h"
#include <assert.h>

int main() {
    auto *p = new Param();
    p->set_alpha(1);
    assert(p->get_alpha()==1);
    p->set_beta(2);
    assert(p->get_beta()==2);
    delete(p);
}