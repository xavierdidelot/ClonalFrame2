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

double Node::get_dist_left() {
    return dist_left;
}

void Node::set_dist_left(double dist_left) {
    this->dist_left = dist_left;
}

double Node::get_dist_right() {
    return dist_right;
}

void Node::set_dist_right(double dist_right) {
    this->dist_right = dist_right;
}

double Node::get_dist_father() {
    return dist_father;
}

Node::Node() {
    left = NULL;
    right = NULL;
    father = NULL;
    dist_father = 0;
    dist_left = 0;
    dist_father = 0;
    name = "";
}

Node::~Node() {
}

string Node::to_string() {
    if (left == NULL) return name;
    return("("+ left ->to_string()+":"+std::to_string(dist_left )+","+
                right->to_string()+":"+std::to_string(dist_right)+")");
}

void Node::set_dist_father(double dist_father) {
    this->dist_father = dist_father;
}
