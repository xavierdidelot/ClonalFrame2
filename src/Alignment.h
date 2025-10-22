//
// Created by Didelot, Xavier on 19/09/2025.
//

#ifndef CF2_ALIGNMENT_H
#define CF2_ALIGNMENT_H

#include <string>
#include <fstream>
#include <vector>
#include "Tree.h"

/**
 * @class Alignment
 * @brief Class representing an alignment
 */
class Alignment {
public:
    ///Destructor of alignment
    ~Alignment();
    ///Constructor of alignment from a fasta file
    explicit Alignment(const std::string &filename);
    ///Constructor of non-polymorphic alignment of n genomes of length l
    Alignment(size_t num,size_t len);
    ///Constructor of alignment from a tree
    explicit Alignment(Tree * t);
    ///Return data for genome i position j
    char get_data(size_t i, size_t j);
    ///Return number of genomes in alignment
    size_t get_n();
    ///Return length of the alignment
    size_t get_l();
    ///Write alignment into a fasta file
    void write(const std::string &filename);
    ///Return name of the i-th genome
    std::string get_name(size_t i);
    ///Pairwise distance
    double dist(size_t i, size_t j);
    ///Return a sequence
    std::vector<uint8_t> * get_seq(size_t i);
    ///Make some of the data be missed at random
    void miss_at_random(double x);

    ///Indicate if a nucleotide is missing
    static inline bool nuc_is_missing(const uint8_t c) {return (c&(1u<<4))!=0;}
    ///Make a nuc from an int (0->A, 1->C, 2->G, 3->T)
    static inline uint8_t nuc_from_int(const int i) {return 1u<<i;}
    ///Any nucleotide
    static inline uint8_t nuc_any() {return 0b00001111;}
    ///Check if a nuc is simply A, C, G or T
    static inline bool nuc_is_single(const uint8_t c) {return c==1u || c==1u<<1 || c==1u<<2 || c==1u<<3;}
    ///Union of nucs
    static inline uint8_t nuc_union(uint8_t a,uint8_t b) {return a|b;}
    ///Pick one if more than one
    static inline uint8_t nuc_pick_one(uint8_t a) {return a & (0 - a);}
    ///Intersect of nucs
    static inline uint8_t nuc_intersect(uint8_t a,uint8_t b) {return a&b;}
    ///Sets the nuc a to nuc b while retaining extra info
    static inline void nuc_set_as(uint8_t * a,const uint8_t b) {*a=((*a)&0b11110000) | (b&0b00001111);}
    ///Indicate if two nucs are the same
    static inline bool nuc_is_same(uint8_t a,uint8_t b) {return (a&0b00001111)==(b&0b00001111);};

private:
    static uint8_t encode(char in);
    static char decode(uint8_t in);
    size_t n;
    size_t l;
    std::vector<std::vector<uint8_t>*> data;
    std::vector<std::string> names;
};


#endif //CF2_ALIGNMENT_H