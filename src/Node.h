//
// Created by Didelot, Xavier on 29/09/2025.
//

#ifndef CLONALFRAME2_NODE_H
#define CLONALFRAME2_NODE_H

#include <string>
#include <iomanip>
#include "global.h"

using std::string,std::vector;

/**
 * @class Node
 * @brief Class representing a node with (potentially) two children nodes
 */
class Node {
public:
    ///Return name of the node
    string get_name() {return name;}
    ///Set name of the node
    void set_name(const string & s) {name = s;}
    ///Return node on the left
    Node * get_left() {return left;}
    ///Set node on the left
    void set_left(Node *l) {left = l;}
    ///Return node on the right
    Node * get_right() {return right;}
    ///Set node on the right
    void set_right(Node * r) {right = r;}
    ///Return father node
    Node * get_father() {return father;}
    ///Set father node
    void set_father(Node * f) {father=f;};
    ///Return distance to father
    double get_dist_father() {return dist_father;}
    ///Set distance to father
    void set_dist_father(double d) {dist_father=d;};
    ///Return sequence
    vector<uint8_t> * get_seq() {return seq;};
    ///Sets sequence
    void set_seq(vector<uint8_t> * s) {seq = s;};
    ///Fitch algorithm, first pass, returns score
    int fitch_pass1(size_t pos);//recursive
    ///Fitch algorithm, second pass
    void fitch_pass2(size_t pos);//recursive

    ///Constructor for a new node
    Node();
    ///Constructor for a new node from a Newick string
    Node(string s,vector <Node*> * nodes);//recursive
    ///Destructor for a node
    ~Node();
    ///Return a newick string for this node and its children
    string to_string();//recursive
    size_t get_index() {return index;};
    void set_index(size_t i) {index = i;};

private:
    string name;
    vector<uint8_t> *seq;
    Node *left;
    Node *right;
    Node *father;
    double dist_father;
    size_t index;
};

#endif //CLONALFRAME2_NODE_H