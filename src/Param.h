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
    double get_theta() {return theta;}
    void set_theta(double theta) {this->theta = theta;}
    double get_rho() {return rho;}
    void set_rho(double rho) {this->rho = rho;}
    double get_delta() {return delta;}
    void set_delta(double delta) {this->delta = delta;}
    double get_nu() {return nu;}
    void set_nu(double nu) {this->nu = nu;}
    ~Param();
    /**
     * Run the main MCMC
     * @param burn Number of iterations to do before sampling
     * @param nbit Total number of iterations to do
     * @param thin Interval between sampling
     * @param outfile File to which output is written
     */
    void mcmc(int burn,int nbit,int thin,string outfile);

private:
    Tree * tree;
    double theta;
    double rho;
    double delta;
    double nu;
};


#endif //CLONALFRAME2_PARAM_H