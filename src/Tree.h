//
// Created by Didelot, Xavier on 29/09/2025.
//

#ifndef CLONALFRAME2_TREE_H
#define CLONALFRAME2_TREE_H

#include <string>
#include <fstream>
#include <vector>
#include "Node.h"

using std::string,std::vector;

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
    static Tree* nj(vector<vector<double>> * d,vector<string> * names);
    ///Create a tree from a Newick string
    static Tree* read_newick(const string &filename);
    ///Create a tree by adding noise to a tree
    static Tree* noise_tree(Tree * t, double eps);
    ///Create a copy of a tree
    static Tree* copy_tree(Tree * t);
    ///Return the tree as a newick string
    string to_string();
    ///Return the tree as a newick string, revealing internal formating
    string to_string_internal();
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
    ///Add mutations on a tree to generate sequence of length l
    int simulate_sequences(size_t l);
    ///Add recombination events and mutations on a tree to generate sequence of length l
    int simulate_sequences_rec(size_t l, double rho, double delta, double nu);
    ///Perform pattern analysis on alignment without ancestors while keeping ancestor sequences
    void simulate_pattern_analysis();
    ///Assign sequences from alignment to the tree leaves
    void set_sequences(Alignment * ali);
    ///Parsimonious ancestral state reconstruction, returning number of substitutions in patterns
    int parsimony();
    ///Parsimonious ancestral state reconstruction, returning number of substitutions in patterns, recursive version
    int parsimony_recursive();
    ///Joint ML ancestral state reconstruction, using Pupko et al algorithm, return likelihood
    double jointml();
    ///Likelihood computation using pruning algorithm
    double likelihood();
    ///Count number of substitutions
    int count_subs();
    ///Sum of branch lengths
    double sum_bra_len();
    ///Likelihood given ancestral sequences
    double likelihood_anc();
    ///Likelihood given ancestral sequences and recombination parameters
    double likelihood_anc_rec(double rho,double delta,double nu);
    ///Reorder nodes
    void reorder();
    ///Compute pairwise distance matrix
    vector<vector<double>> * distmat();
    ///Return alignment
    Alignment * get_alignment()  {return ali;}
    ///Perform SPR move and return log of the proposal probability
    double rSPR(string * name_moving);
    ///Select new bro for SPR move
    double select_rSPR(size_t moving, size_t * newbro, double * m, double * f, double * b, double * d);
    ///Tree prior
    double prior();
    ///Probability of proposing from this to tree, given name of moving node
    double lq_rSPR(Tree *tree2, string name_moving);
    ///dgamma function, returning log(pdf)
    static inline double dgamma(double x,double shape,double scale)
        {return (shape - 1.0) * log(x) - x / scale - shape * log(scale) - lgamma(shape);}
    ///Determine the shape and scale parameters of the gamma distribution for the length of a branch
    void determine_shape_scale(double diffs, double * shape, double * scale);
    ///Count number of differences between ancestral sequences of two nodes
    int count_diffs(Node * n1, Node * n2);

private:
    Tree() = default;
    size_t n;
    vector <Node*> nodes;
    Alignment * ali;
};


#endif //CLONALFRAME2_TREE_H