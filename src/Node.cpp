//
// Created by Didelot, Xavier on 29/09/2025.
//

#include "Node.h"


string Node::get_name() {
    return name;
}

void Node::set_name(string name) {
    this->name = name;
}

Node * Node::get_left() {
    return left;
}

void Node::set_left(Node *left) {
    this->left = left;
}

Node * Node::get_right() {
    return right;
}

void Node::set_right(Node *right) {
    this->right = right;
}

Node * Node::get_father() {
    return father;
}

void Node::set_father(Node *father) {
    this->father = father;
}

double Node::get_dist_father() {
    return dist_father;
}

Node::Node() {
    left = NULL;
    right = NULL;
    father = NULL;
    dist_father = 0;
    name = "";
    seq=NULL;
}

Node::~Node() {
}

string Node::to_string() {
    if (left == NULL) return name;
    return("("+ left ->to_string()+":"+std::to_string(left ->get_dist_father())+","+
                right->to_string()+":"+std::to_string(right->get_dist_father())+")");
}

void Node::add_mutations(string * above, double theta) {
    seq=new string(*above);
    int l=seq->size();
    uniform_int_distribution unif4(0,3);
    uniform_int_distribution unifL(0,l-1);
    poisson_distribution pois(dist_father*theta/2);
    int m=pois(rng);
    for (int i=0;i<m;i++) {
        int p=unifL(rng);
        while (above->at(p)==seq->at(p)) {
            int w=unif4(rng);
            *(seq->data()+p)=*(lets+w);
        }
    }
    if (left!= NULL) left ->add_mutations(seq,theta);
    if (right!=NULL) right->add_mutations(seq,theta);
}

void Node::set_dist_father(double dist_father) {
    this->dist_father = dist_father;
}

string * Node::get_seq() {
    return seq;
}
