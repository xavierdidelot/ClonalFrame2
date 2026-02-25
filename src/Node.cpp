//
// Created by Didelot, Xavier on 29/09/2025.
//

#include "Node.h"
#include "Alignment.h"

using std::string,std::cerr,std::cout,std::endl,std::vector;

string Node::get_name() {
    return name;
}

void Node::set_name(const string &s) {
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

Node::Node(const string &s, vector <Node*> * nodes) {
    left=nullptr;
    right=nullptr;
    father=nullptr;
    seq=nullptr;
    dist_father=0;
    if (s.front()=='(') {
        size_t last=s.size()-1;
        while (s[last]!=')') last--;
        if (last<s.size()-1) {
            dist_father=stof(s.substr(last+2));
        }
        size_t mid=1;
        int counter=0;
        while (true) {
            if (s[mid]=='(') counter++;
            if (s[mid]==')') counter--;
            if (s[mid]==',' && counter==0) break;
            mid++;
        }
        nodes->push_back(this);
        string s1=s.substr(1,mid-1);
        string s2=s.substr(mid+1,last-mid-1);
        left=new Node(s1, nodes);
        right=new Node(s2, nodes);
        left->set_father(this);
        right->set_father(this);
    } else {
        nodes->insert(nodes->begin(),this);
        std::istringstream iss(s);
        string tmp;
        getline(iss,tmp,':');
        set_name(tmp);
        getline(iss,tmp);
        set_dist_father(stof(tmp));
    }
}

Node::~Node() {
    delete seq;
}

string Node::to_string() {
    if (left == nullptr) return name;
    std::ostringstream oss;
    oss << "(";
    oss << left->to_string();
    oss << std::scientific << std::setprecision(6) << ":" << left ->get_dist_father();
    oss << ",";
    oss << right->to_string();
    oss << std::scientific << std::setprecision(6) << ":" << right->get_dist_father();
    oss << ")";
    return oss.str();
}

void Node::simulate_sequences(vector<uint8_t> * above) {
    seq=new vector<uint8_t>(*above);
    int l=seq->size();
    std::uniform_int_distribution unif4(0,3);
    std::uniform_int_distribution unif_l(0,l-1);
    std::poisson_distribution pois(dist_father*l);
    int m=pois(rng);
    for (int i=0;i<m;i++) {
        int p=unif_l(rng);
        while (above->at(p)==seq->at(p)) {
            seq->at(p)=Alignment::nuc_from_int(unif4(rng));
        }
    }
    if (left!= nullptr) left ->simulate_sequences(seq);
    if (right!=nullptr) right->simulate_sequences(seq);
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
        if (Alignment::nuc_is_unsure(seq->at(pos))) seq->at(pos)=0b00011111;
        return 0;
    }
    int s=left->fitch_pass1(pos);
    s+=right->fitch_pass1(pos);
    uint8_t toset=left->seq->at(pos)&right->seq->at(pos)&0b00001111;//intersect
    if (toset==0u) {
        s++;
        toset=left->seq->at(pos)|right->seq->at(pos);//union
    }
    Alignment::nuc_set_as(&seq->at(pos),toset);
    return s;
}

void Node::fitch_pass2(size_t pos) {
    uint8_t toset;
    if (father==nullptr) toset=seq->at(pos)&0b00001111;
    else {
        toset=seq->at(pos)&father->seq->at(pos)&0b00001111;
        if (toset==0u) toset=seq->at(pos)&0b00001111;
    }
    toset=toset&-toset;//pick first common bit
    Alignment::nuc_set_as(&seq->at(pos),toset);//preserve missing status
    if (left!=nullptr) {
        left ->fitch_pass2(pos);
        right->fitch_pass2(pos);
    }
}