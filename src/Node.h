//
// Created by Didelot, Xavier on 29/09/2025.
//

#ifndef CLONALFRAME2_NODE_H
#define CLONALFRAME2_NODE_H

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "RNGGlobal.h"

using namespace std;

/**
 * @class Node
 * @brief Class representing a node with (potentially) two children nodes
 */
class Node {
public:
    ///Return name of the node
    string get_name();
    ///Set name of the node
    void set_name(const string &);
    ///Return node on the left
    Node * get_left();
    ///Set node on the left
    void set_left(Node *);
    ///Return node on the right
    Node * get_right();
    ///Set node on the right
    void set_right(Node *);
    ///Return father node
    Node * get_father();
    ///Set father node
    void set_father(Node *);
    ///Return distance to father
    double get_dist_father();
    ///Set distance to father
    void set_dist_father(double);
    ///Return sequence
    vector<uint8_t> * get_seq();
    ///Sets sequence
    void set_seq(vector<uint8_t> *);
    ///Fitch algorithm, first pass, returns score
    int fitch_pass1(size_t pos);
    ///Fitch algorithm, second pass
    void fitch_pass2(size_t pos);

    ///Constructor for a new node
    Node();
    ///Constructor for a new node from a Newick string
    Node(const string &s);
    ///Destructor for a node
    ~Node();
    ///Return a newick string for this node and its children
    string to_string();
    ///Add mutations
    void add_mutations(vector<uint8_t> * above,double theta);

private:
    string name;
    vector<uint8_t> *seq;
    Node *left;
    Node *right;
    Node *father;
    double dist_father;
};


#endif //CLONALFRAME2_NODE_H