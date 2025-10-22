//
// Created by Didelot, Xavier on 21/10/2025.
//

#include "Param.h"

using std::string,std::cerr,std::cout,std::endl,std::vector;

Param::Param(Tree *t) {
    tree = t;
    theta=10;
    rho=theta/10;
    nu=0.01;
    delta=0.001;
}

Param::~Param() {
}

void Param::mcmc(int burn, int nbit, int thin) {
    cout<<burn<<" "<<nbit<<" "<<thin<<endl;
}
