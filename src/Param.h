//
// Created by Didelot, Xavier on 21/10/2025.

#ifndef CLONALFRAME2_PARAM_H
#define CLONALFRAME2_PARAM_H

#include "Tree.h"

/**
 * Class containing all parameters
 */
class Param {
public:
    ///Constructor given a tree
    explicit Param(Tree * t);
    ~Param();
    /**
     * Run the main MCMC
     * @param burn Number of iterations to do before sampling
     * @param nbit Total number of iterations to do
     * @param thin Interval between sampling
     */
    void mcmc(int burn,int nbit,int thin);

private:
    Tree * tree;
    double theta;
    double rho;
    double delta;
    double nu;
};


#endif //CLONALFRAME2_PARAM_H