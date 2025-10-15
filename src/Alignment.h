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
    Alignment(int num,int len);
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
    ///Return a sequence
    vector<uint8_t> * get_seq(int i);
    ///Indicate if a nucleotide is missing
    static inline bool nuc_is_missing(const uint8_t c) {return (c&(1u<<4))!=0;}
    static inline uint8_t nuc_from_int(const int i) {return 1u<<i;}
    static inline bool nuc_is_single(const uint8_t c) {return c==1u || c==1u<<1 || c==1u<<2 || c==1u<<3;}
    static inline uint8_t nuc_or(uint8_t a,uint8_t b) {return a|b;}
    static inline void nuc_set_as(uint8_t * a,const uint8_t b) {*a=((*a)&0b11110000) | (b&0b00001111);}

private:
    static uint8_t encode(char in);
    static char decode(uint8_t in);
    int n;
    int l;
    vector<vector<uint8_t>*>* data;
    vector<string> names;
};


#endif //CF2_ALIGNMENT_H