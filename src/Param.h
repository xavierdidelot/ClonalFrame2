//
// Created by Didelot, Xavier on 21/10/2025.
//

#ifndef CLONALFRAME2_PARAM_H
#define CLONALFRAME2_PARAM_H

#include "Tree.h"

class Param {
public:
    Param(Tree * t);
    ~Param() = default;
    void mcmc(int burn,int nbit,int thin);
private:
    Tree * tree;
    double theta;
    double rho;
    double delta;
    double nu;
};


#endif //CLONALFRAME2_PARAM_H