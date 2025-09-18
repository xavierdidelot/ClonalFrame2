//
// Created by Didelot, Xavier on 17/09/2025.
//

#ifndef UNTITLED_PARAM_H
#define UNTITLED_PARAM_H
#include <iostream>
#include <ostream>


class Param {
private:
    double alpha;
    double beta;
    double gamma;

public:
    double get_alpha() const;

    double get_beta() const;

    double get_gamma() const;

    void set_alpha(double alpha);

    void set_beta(double beta);

    void set_gamma(double gamma);

    void print();

    Param(double alpha, double beta, double gamma);

    Param();
};


#endif //UNTITLED_PARAM_H