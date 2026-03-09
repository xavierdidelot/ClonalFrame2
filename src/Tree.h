//
// Created by Didelot, Xavier on 29/09/2025.
//

#ifndef CLONALFRAME2_TREE_H
#define CLONALFRAME2_TREE_H

#include <string>
#include <fstream>
#include <vector>
#include "Node.h"

using std::string,std::cerr,std::cout,std::endl,std::vector;

class Alignment;

/**
 * @class Tree
 * @brief Class representing an unrooted non-ultrametric tree
 */
class Tree {
public:
    ///Create a deterministic tree with n leaves
    static Tree* deterministic_tree(size_t n);
    ///Create a random tree with n leaves
    static Tree* random_tree(size_t n,double avg_bra_len);
    ///Create a tree using the NJ algorithm
    static Tree* nj(Alignment * ali);
    ///Create a tree from a Newick string
    static Tree* read_newick(const string &filename);
    ///Return the tree as a newick string
    string to_string();
    ///Destructor for a tree
    ~Tree();
    ///Return number of leaves in the tree
    size_t get_n() {return n;}
    ///Return node with a given index
    Node * get_node(size_t i) {return nodes[i];}
    ///Write tree as a newick file
    void write(const string &filename);
    ///Check that the tree structure is correct, especially forward/backward links
    bool check_structure();
    ///Add mutations on a tree to generate sequence of length l with mutation rate theta/2
    int simulate_sequences(size_t l);
    ///Assign sequences from alignment to the tree leaves
    void set_sequences(Alignment * ali);
    ///Parsimonious ancestral state reconstruction, returning number of substitutions in patterns
    int parsimony();
    ///Parsimonious ancestral state reconstruction, returning number of substitutions in patterns, recursive version
    double jointml();
    ///Joint ML ancestral state reconstruction, using Pupko et al algorithm, return likelihood
    int parsimony_recursive();
    ///Count number of substitutions in patterns
    int count_subs();
    ///Sum of branch lengths
    double sum_bra_len();
    ///Likelihood given ancestral sequences
    double likelihood_anc();
    ///Reorder nodes
    void reorder();
    size_t get_npat() {return npat;}
    vector<size_t> * get_patindex()  {return patindex;}
    vector<size_t> * get_weights()  {return weights;}

private:
    Tree() = default;
    size_t n;
    vector <Node*> nodes;
    size_t npat;
    vector<size_t>* patindex;
    vector<size_t>* weights;
};


#endif //CLONALFRAME2_TREE_H