//
// Created by Didelot, Xavier on 21/10/2025.

#ifndef CLONALFRAME2_PARAM_H
#define CLONALFRAME2_PARAM_H

#include "Tree.h"
#include <iostream>

/**
 * Class containing all parameters
 */
class Param {
public:
    ///Constructor given a tree
    explicit Param(Tree * t);
    double get_rho() {return rho;}
    void set_rho(double rho) {this->rho = rho;}
    double get_delta() {return delta;}
    void set_delta(double delta) {this->delta = delta;}
    double get_nu() {return nu;}
    void set_nu(double nu) {this->nu = nu;}
    Tree * get_tree() {return tree;}
    ~Param();
    /**
     * Run the main MCMC
     * @param burn Number of iterations to do before sampling
     * @param nbit Total number of iterations to do
     * @param thin Interval between sampling
     * @param outfile File to which output is written
     */
    void mcmc(int burn,int nbit,int thin);
    ///Write to json file
    void to_json_file();
    ///Selection of updates to perform in MCMC
    void select_updates(const string & string);

private:
    void print_progress(int current, int total);
    Tree * tree;
    double rho;
    double delta;
    double nu;
    double ll;
    double tree_prior;
    vector<double> paccepted;
    int update_tree;
    int update_rho;
    int update_delta;
    int update_nu;
    int update_lengths;
};


#endif //CLONALFRAME2_PARAM_H