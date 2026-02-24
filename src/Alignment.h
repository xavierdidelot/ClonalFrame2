//
// Created by Didelot, Xavier on 19/09/2025.
//

#ifndef CF2_ALIGNMENT_H
#define CF2_ALIGNMENT_H

#include <string>
#include <fstream>
#include <vector>
#include "Tree.h"

using std::string,std::cerr,std::cout,std::endl,std::vector;

/*Bit encoding of nucleotide (from right to left)
Bit1=is_currently(A)
Bit2=is_currently(C)
Bit3=is_currently(G)
Bit4=is_currently(T)
Bit5=is_unsure
Bit6=unused
Bit7=unused
Bit8=unused
*/

/**
 * @class Alignment
 * @brief Class representing an alignment
 */
class Alignment {
public:
    ///Create an alignment from a fasta file
    static Alignment* read_fasta(const string &filename);
    ///Constructor of non-polymorphic alignment of n genomes of length l
    Alignment(size_t num,size_t len);
    ///Constructor of alignment from a tree
    explicit Alignment(Tree * t);
    ///Destructor of alignment
    ~Alignment();
    ///Return data for genome i position j
    char get_data(size_t i, size_t j);
    ///Return number of genomes in alignment
    size_t get_n();
    ///Return length of the alignment
    size_t get_l();
    ///Return number of patterns
    size_t get_npat();
    ///Write alignment into a fasta file
    void write(const string &filename);
    ///Return name of the i-th genome
    string get_name(size_t i);
    ///Pairwise distance
    double dist(size_t i, size_t j);
    ///Return a sequence
    vector<uint8_t> * get_seq(size_t i);
    ///Make some of the data be missed at random
    void miss_at_random(double x);
    ///Pattern analysis
    void pattern_analysis();

    ///Indicate if a nucleotide is unsure
    static inline bool nuc_is_unsure(const uint8_t c) {return (c&0b00010000)!=0;}
    ///Make a nuc from an int (0=A, 1=C, 2=G, 3=T, 4=N)
    static inline uint8_t nuc_from_int(const int i) {return 1u<<i;}
    ///Sets the nuc a to nuc b while retaining extra info
    static inline void nuc_set_as(uint8_t * a,const uint8_t b) {*a=((*a)&0b11110000) | (b&0b00001111);}
    ///Indicate if two nucs are currently the same
    static inline bool nuc_is_currently_same(uint8_t a,uint8_t b) {return (a&0b00001111)==(b&0b00001111);};

private:
    Alignment() = default;
    static uint8_t encode(char in);
    static char decode(uint8_t in);
    size_t n;
    size_t l;
    vector<vector<uint8_t>*> data;
    vector<string> names;
    size_t npat;
    vector<size_t> patindex;
    vector<size_t> weights;
};


#endif //CF2_ALIGNMENT_H