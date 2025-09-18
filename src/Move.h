//
// Created by Didelot, Xavier on 17/09/2025.
//

#ifndef UNTITLED_MOVE_H
#define UNTITLED_MOVE_H
#include "Param.h"


class Move {
private:
    double byhowmuch;
public:
    Move(double byhowmuch);

    void action(Param *p);
};


#endif //UNTITLED_MOVE_H