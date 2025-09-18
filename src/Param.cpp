//
// Created by Didelot, Xavier on 17/09/2025.
//

#include "Param.h"

double Param::get_alpha() const {
    return alpha;
}

double Param::get_beta() const {
    return beta;
}

double Param::get_gamma() const {
    return gamma;
}

void Param::set_alpha(double alpha) {
    this->alpha = alpha;
}

Param::Param(double alpha, double beta, double gamma) {
    this->alpha = alpha;
    this->beta = beta;
    this->gamma = gamma;
}

Param::Param() {}

void Param::set_beta(double beta) {
    this->beta = beta;
}

void Param::set_gamma(double gamma) {
    this->gamma = gamma;
}

void Param::print() {
    std::cout<<"alpha="<<alpha<<", beta="<<beta<<", gamma="<<gamma<<std::endl;
}
