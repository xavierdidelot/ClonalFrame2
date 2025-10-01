//
// Created by Didelot, Xavier on 29/09/2025.
//

#include "Tree.h"
#include "Alignment.h"

#include <iostream>

Tree::Tree() {
}

Tree::Tree(int n) {
    this->n=n;
    nodes=vector<Node*>();
    for (int i = 0; i < n; i++) {
        Node * node = new Node();
        node->set_name("seq"+std::to_string(i+1));
        node->set_dist_father(1.0);
        nodes.push_back(node);
    }
    from_root = vector<Node*>();
    from_root.push_back(nodes[0]);
    from_root.push_back(nodes[1]);
    from_root.push_back(nodes[n-1]);

    for (int i = n-2; i >= 2; i--) {
        Node * node = new Node();
        nodes.push_back(node);
        node->set_left(nodes[i]);
        nodes[i]->set_father(node);
        node->set_dist_father(1.0);
        node->set_right(from_root[2]);
        from_root[2]->set_father(node);
        from_root[2]=node;
    }

}

Tree::Tree(int n, double avg_bra_len) {
    this->n=n;
    nodes=vector<Node*>();
    vector<Node*> to_merge;
    to_merge=vector<Node*>();
    exponential_distribution rexp(avg_bra_len);
    for (int i = 0; i < n; i++) {
        Node * node = new Node();
        node->set_name("seq"+std::to_string(i+1));
        node->set_dist_father(rexp(rng));
        nodes.push_back(node);
        to_merge.push_back(node);
    }
    from_root = vector<Node*>();
    while (to_merge.size() > 3) {
        Node * node = new Node();
        nodes.push_back(node);
        uniform_int_distribution unif(0, (int)(to_merge.size()-1));
        int a=unif(rng);
        int b=a;
        while (a==b) b=unif(rng);
        node->set_left(to_merge[a]);
        node->set_right(to_merge[b]);
        to_merge[a]->set_father(node);
        to_merge[b]->set_father(node);
        node->set_dist_father(rexp(rng));
        to_merge[a]=node;
        to_merge[b]=to_merge.back();
        to_merge.pop_back();
    }
    for (int i=0;i<3;i++) from_root.push_back(to_merge[i]);
}

Tree::Tree(Alignment *a) {
    int n=a->get_n();
    this->n=n;
    nodes=vector<Node*>();
    for (int i = 0; i < n; i++) {
        Node * node = new Node();
        node->set_name(a->get_name(i));
    }
    //TODO
}

string Tree::to_string() {
    return("("+ from_root[0]->to_string()+":"+std::to_string(from_root[0]->get_dist_father())+","+
                from_root[1]->to_string()+":"+std::to_string(from_root[1]->get_dist_father())+","+
                from_root[2]->to_string()+":"+std::to_string(from_root[2]->get_dist_father())+");");
}

Tree::~Tree() {
}

int Tree::get_n() {
    return n;
}

Node * Tree::get_node(int i) {
    return nodes[i];
}

void Tree::write(string filename) {
    ofstream file;
    file.open(filename);
    if (file.is_open()==false) {
        cout<<"Error opening file "<<filename<<endl;
        exit(1);
    }
    file<<to_string()<<endl;
}

bool Tree::check_structure() {
    if (this->n*2-3 - nodes.size()!=0) return false;
    if (from_root.size()!=3) return false;
    for (int i=0;i<3;i++) if (from_root[i]->get_father()!=NULL) return false;
    for (int i=0; i<nodes.size(); i++) {
        if (nodes[i]->get_left ()!=NULL && nodes[i]->get_left ()->get_father()!=nodes[i]) return false;
        if (nodes[i]->get_right()!=NULL && nodes[i]->get_right()->get_father()!=nodes[i]) return false;
    }
    return true;
}

void Tree::add_mutations(int l, double theta) {
    seq=new string(l,'A');
    uniform_int_distribution unif4(0,3);
    for (int i=0;i<l;i++) {
        *(seq->data()+i)=*(lets+unif4(rng));
    }
    for (int i=0;i<3;i++) from_root[i]->add_mutations(seq,theta);

}
