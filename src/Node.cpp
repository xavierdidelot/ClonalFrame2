//
// Created by Didelot, Xavier on 29/09/2025.
//

#include "Node.h"
#include "Alignment.h"

Node::Node() {
    left = nullptr;
    right = nullptr;
    father = nullptr;
    dist_father = 0;
    seq=nullptr;
}

Node::Node(string s, vector <Node*> * nodes) {
    left=nullptr;
    right=nullptr;
    father=nullptr;
    seq=nullptr;
    dist_father=0;
    if (s.front()=='(') {
        //Find last closing bracket
        size_t last=s.size()-1;
        while (s[last]!=')') last--;
        //Test for multifurcation and if so add brackets
        vector<int> commas=vector<int>();
        int counter=0;
        for (size_t cur=1;cur<s.size()-1;cur++) {
            if (s[cur]=='(') counter++;
            if (s[cur]==')') counter--;
            if (s[cur]==',' && counter==0) commas.push_back(cur);
        }
        if (commas.size()>1) {
            s.insert(last-1,"):0.0");
            s.insert(commas[0]+1,"(");
        }

        std::istringstream iss(s.substr(last+1));
        string b;
        getline(iss,b,':');
        set_name(b);
        if (iss.peek() != EOF) {
            getline(iss,b);
            set_dist_father(stof(b));
        }
        size_t mid=1;
        counter=0;
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
        string b;
        getline(iss,b,':');
        set_name(b);
        getline(iss,b);
        set_dist_father(stof(b));
    }
}

string Node::to_string() {
    std::ostringstream oss;
    if (left == nullptr) std::print(oss,"{}:{:.10e}",name,get_dist_father());
    else std::print(oss,"({},{}){}:{:.10e}",left->to_string(),right->to_string(),name,get_dist_father());
    return oss.str();
}

int Node::simulate_sequences(vector<uint8_t> * above) {
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
    if (left!= nullptr) m+=left ->simulate_sequences(seq);
    if (right!=nullptr) m+=right->simulate_sequences(seq);
    return m;
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
    Alignment::nuc_set_as(&seq->at(pos),toset);//preserve missing status
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