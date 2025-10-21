//
// Created by Didelot, Xavier on 21/10/2025.
//

#include "Param.h"

Param::Param(Tree *t) {
    tree = t;
    theta=10;
    rho=theta/10;
    nu=0.01;
    delta=0.001;
}

void Param::mcmc(int burn, int nbit, int thin) {
    swap(burn,burn);
    swap(nbit,nbit);
    swap(thin,thin);
}
