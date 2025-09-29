//
// Created by Didelot, Xavier on 29/09/2025.
//

#include "Tree.h"

#include <iostream>

Tree::Tree() {
}

Tree::Tree(int n) {
    this->n=n;
    for (int i = 0; i < n; i++) {
        Node * node = new Node();
        node->set_name("seq"+std::to_string(i));
        node->set_dist_father(1.0);
        nodes.push_back(node);
    }
    from_root_1 = nodes[0];
    from_root_2 = nodes[1];
    from_root_3 = nodes[2];

    for (int i = 3; i < n; i++) {
        Node * node = new Node();
        nodes.push_back(node);
        node->set_left(nodes[i]);
        node->set_dist_father(1.0);
        node->set_dist_left(1.0);
        node->set_dist_right(1.0);
        node->set_right(from_root_3);
        from_root_3=node;
    }

}

string Tree::to_string() {
    return("("+ from_root_1->to_string()+":"+std::to_string(from_root_1->get_dist_father())+","+
                from_root_2->to_string()+":"+std::to_string(from_root_2->get_dist_father())+","+
                from_root_3->to_string()+":"+std::to_string(from_root_3->get_dist_father())+")");
}

Tree::~Tree() {
}

int Tree::get_n() {
    return n;
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
