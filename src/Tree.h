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

class Tree {
public:
    Tree();
    Tree(int n);
    string to_string();
    ~Tree();
    int get_n();
    void write(string filename);
private:
    int n;
    Node* from_root_1;
    Node* from_root_2;
    Node* from_root_3;
    vector <Node*> nodes;
};


#endif //CLONALFRAME2_TREE_H