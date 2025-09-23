//
// Created by Didelot, Xavier on 19/09/2025.
//

#ifndef CF2_ALIGNMENT_H
#define CF2_ALIGNMENT_H

#include <string>
#include <fstream>
#include <gsl/gsl_matrix.h>
using namespace std;

class Alignment {
public:
    ~Alignment();

    Alignment(string filename);

    char standardise(char in);

    char getData(int i, int j);

    int getN();

    int getL();

private:
    gsl_matrix_char* data;
};


#endif //CF2_ALIGNMENT_H