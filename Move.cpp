//
// Created by Didelot, Xavier on 17/09/2025.
//

#include "Move.h"

Move::Move(double byhowmuch) {
    this->byhowmuch = byhowmuch;
}

void Move::action(Param *p) {
    p->set_alpha(p->get_alpha()+byhowmuch);
}
