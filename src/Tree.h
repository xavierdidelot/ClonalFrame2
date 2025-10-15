//
// Created by Didelot, Xavier on 29/09/2025.
//

#ifndef CLONALFRAME2_TREE_H
#define CLONALFRAME2_TREE_H

#include <string>
#include <fstream>
#include <vector>
#include "Node.h"

using namespace std;

class Alignment;

/**
 * @class Tree
 * @brief Class representing an unrooted non-ultrametric tree
 */
class Tree {
public:
    ///Constructor for a deterministic tree with n leaves
    explicit Tree(int n);
    ///Return the tree as a newick string
    string to_string();
    ///Destructor for a tree
    ~Tree();
    ///Return number of leaves in the tree
    int get_n();
    Node * get_node(int i);
    ///Write tree as a newick file
    void write(string filename);
    ///Check that the tree structure is correct, especially forward/backward links
    bool check_structure();
    ///Add mutations on a tree to generate sequence of length l with mutation rate theta/2
    void add_mutations(int l,double theta);
    ///Randomise the tree
    void randomise(double avg_bra_len);
    ///NJ algorithm
    void nj(Alignment * ali);
    ///Assign sequences from alignment to the tree leaves
    void set_seqs(Alignment * ali);
    ///Ancestral state reconstruction
    void ancestral_reconstruction();

private:
    int n;
    vector <Node*> nodes;
};


#endif //CLONALFRAME2_TREE_H