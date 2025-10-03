//
// Created by Didelot, Xavier on 19/09/2025.
//

#ifndef CF2_ALIGNMENT_H
#define CF2_ALIGNMENT_H

#include <string>
#include <fstream>
#include <vector>
#include "Tree.h"
using namespace std;

/**
 * @class Alignment
 * @brief Class representing an alignment
 */
class Alignment {
public:
    ///Destructor of alignment
    ~Alignment();
    ///Constructor of alignment from a fasta file
    explicit Alignment(string filename);
    ///Constructor of non-polymorphic alignment of n genomes of length l
    Alignment(int n,int l);
    ///Constructor of alignment from a tree
    explicit Alignment(Tree * t);
    ///Return data for genome i position j
    char get_data(int i, int j);
    ///Return number of genomes in alignment
    int get_n();
    ///Return length of the alignment
    int get_l();
    ///Write alignment into a fasta file
    void write(string filename);
    ///Return name of the i-th genome
    string get_name(int i);
    ///Pairwise distance
    double dist(int i, int j);

private:
    char standardise(char in);
    int n;
    int l;
    vector<string> data;
    vector<string> names;
};


#endif //CF2_ALIGNMENT_H