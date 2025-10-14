//
// Created by Didelot, Xavier on 29/09/2025.
//

#ifndef CLONALFRAME2_NODE_H
#define CLONALFRAME2_NODE_H

#include <string>
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
    void set_name(string);
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
    string * get_seq();
    ///Sets sequence
    void set_seq(string *);
    ///Constructor for a new node
    Node();
    ///Destructor for a node
    ~Node();
    ///Return a newick string for this node and its children
    string to_string();
    ///Add mutations
    void add_mutations(string * above,double theta);

private:
    string name;
    string *seq;
    Node *left;
    Node *right;
    Node *father;
    double dist_father;
};


#endif //CLONALFRAME2_NODE_H