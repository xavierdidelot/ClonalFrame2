//
// Created by Didelot, Xavier on 29/09/2025.
//

#ifndef CLONALFRAME2_NODE_H
#define CLONALFRAME2_NODE_H

#include <string>
using namespace std;

class Node {
public:
    string get_name();
    void set_name(string name);
    Node * get_left();
    void set_left(Node *left);
    Node * get_right();
    void set_right(Node *right);
    Node * get_father();
    void set_father(Node *father);
    double get_dist_father();
    void set_dist_father(double dist_father);
    Node();
    ~Node();
    string to_string();

private:
    string name;
    Node *left;
    Node *right;
    Node *father;
    double dist_father;
};


#endif //CLONALFRAME2_NODE_H