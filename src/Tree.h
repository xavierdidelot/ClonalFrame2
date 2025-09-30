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

/**
 * @class Tree
 * @brief Class representing an unrooted non-ultrametric tree
 */
class Tree {
public:
    ///Constructor for an empty tree
    Tree();
    ///Constructor for a tree with n leaves and deterministic structure
    Tree(int n);
    ///Return the tree as a newick string
    string to_string();
    ///Destructor for a tree
    ~Tree();
    ///Return number of leaves in the tree
    int get_n();
    ///Write tree as a newick file
    void write(string filename);
    ///Check that the tree structure is correct, especially forward/backward links
    bool check_structure();
private:
    int n;
    Node* from_root_1;
    Node* from_root_2;
    Node* from_root_3;
    vector <Node*> nodes;
};


#endif //CLONALFRAME2_TREE_H