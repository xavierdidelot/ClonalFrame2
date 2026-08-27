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
     * @param nbit Total number of iterations to do
     * @param thin Interval between sampling
     */
    void mcmc(int nbit,int thin);
    void adapt(int index,double mh,int i,int rep);

    ///Write to json file
    void state_to_json_file(int iter);
    ///Return pointer to vector of updates
    vector<int> * get_updates() {return &updates;}
    ///Return pointer to vector of prior params
    vector<double> * get_prior_params() {return &prior_params;}

private:
    void print_progress(int current, int total);
    Tree * tree;
    double rho;
    double delta;
    double nu;
    double ll;
    double prior_tree;
    double prior_rho;
    double prior_delta;
    double prior_nu;
    vector<int> naccepted;
    vector<int> updates;
    vector<double> sd;
    const double adapt_target=0.234;
    const double adapt_scale=1.0/adapt_target/(1.0-adapt_target);
    vector<double> adapt_nstart;
    vector<double> adapt_lambda;
    vector<double> prior_params;
};


#endif //CLONALFRAME2_PARAM_H