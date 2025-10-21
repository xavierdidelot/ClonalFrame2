//
// Created by Didelot, Xavier on 29/09/2025.
//

#include "Node.h"

#include "Alignment.h"


string Node::get_name() {
    return name;
}

void Node::set_name(string s) {
    name = s;
}

Node * Node::get_left() {
    return left;
}

void Node::set_left(Node *l) {
    left=l;
}

Node * Node::get_right() {
    return right;
}

void Node::set_right(Node *r) {
    right=r;
}

Node * Node::get_father() {
    return father;
}

void Node::set_father(Node *f) {
    father=f;
}

double Node::get_dist_father() {
    return dist_father;
}

Node::Node() {
    left = nullptr;
    right = nullptr;
    father = nullptr;
    dist_father = 0;
    seq=nullptr;
}

string Node::to_string() {
    if (left == nullptr) return name;
    ostringstream oss;
    oss << "(";
    oss << left->to_string();
    oss << scientific << setprecision(6) << ":" << left ->get_dist_father();
    oss << ",";
    oss << right->to_string();
    oss << scientific << setprecision(6) << ":" << right->get_dist_father();
    oss << ")";
    return oss.str();
}

void Node::add_mutations(vector<uint8_t> * above, double theta) {
    seq=new vector<uint8_t>(*above);
    int l=seq->size();
    uniform_int_distribution unif4(0,3);
    uniform_int_distribution unif_l(0,l-1);
    poisson_distribution pois(dist_father*theta/2);
    int m=pois(rng);
    for (int i=0;i<m;i++) {
        int p=unif_l(rng);
        while (above->at(p)==seq->at(p)) {
            seq->at(p)=Alignment::nuc_from_int(unif4(rng));
        }
    }
    if (left!= nullptr) left ->add_mutations(seq,theta);
    if (right!=nullptr) right->add_mutations(seq,theta);
}

void Node::set_dist_father(double d) {
    dist_father=d;
}

vector<uint8_t> * Node::get_seq() {
    return seq;
}

void Node::set_seq(vector<uint8_t> *s) {
    seq=s;
}

int Node::fitch_pass1(size_t pos) {
    if (left==nullptr) {
        if (Alignment::nuc_is_missing(seq->at(pos))) Alignment::nuc_set_as(&(seq->at(pos)),Alignment::nuc_any());
        return 0;
    }
    int s=left->fitch_pass1(pos);
    s+=right->fitch_pass1(pos);
    uint8_t inter=Alignment::nuc_intersect(left->get_seq()->at(pos),right->get_seq()->at(pos));
    if (inter==0u) {
        s++;
        seq->at(pos)=Alignment::nuc_union(left->get_seq()->at(pos),right->get_seq()->at(pos));
    } else {
        seq->at(pos)=inter;
    }
    return s;
}

void Node::fitch_pass2(size_t pos) {
    uint8_t f;
    if (father==nullptr) f=seq->at(pos);
    else {
        f=Alignment::nuc_intersect(get_seq()->at(pos),father->get_seq()->at(pos));
        if (f==0u) f=get_seq()->at(pos);
    }
    seq->at(pos)=Alignment::nuc_pick_one(f);
    if (left!=nullptr) {
        left ->fitch_pass2(pos);
        right->fitch_pass2(pos);
    }
}