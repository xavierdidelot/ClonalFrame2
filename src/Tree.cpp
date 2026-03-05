//
// Created by Didelot, Xavier on 29/09/2025.

#include "Tree.h"
#include "Alignment.h"

#include <iostream>

using std::string,std::cerr,std::cout,std::endl,std::vector;

Tree * Tree::deterministic_tree(size_t n) {
    auto t = new Tree();
    t->n=n;
    t->nodes=vector<Node *>(n+n-1);
    for (size_t i = 0; i < n; i++) {
        auto node = new Node();
        node->set_name("seq" + std::to_string(i + 1));
        node->set_dist_father(1.0);
        t->nodes[i]=node;
    }

    for (size_t i = 0; i < n-1; i++) {
        auto node = new Node();
        t->nodes[n+i]=node;
        node->set_left(t->nodes[i]);
        t->nodes[i]->set_father(t->nodes[n+i]);
        node->set_dist_father(1.0);
        if (i>0) {
            node->set_father(t->nodes[n+i-1]);
            t->nodes[n+i-1]->set_right(t->nodes[n+i]);
        }
    }
    t->nodes[n+n-2]->set_right(t->nodes[n-1]);
    t->nodes[n-1]->set_father(t->nodes[n+n-2]);
    t->nodes[n-1]->set_dist_father(1.0);
    return t;
}

Tree* Tree::read_newick(const string &filename) {
    auto t=new Tree();
    string line;
    std::ifstream file;
    if (verbose) cout << "Opening file " << filename << endl;
    file.open(filename);
    if (file.is_open() == false) {
        cerr << "Error opening file " << filename << endl;
        exit(1);
    }
    getline(file, line);
    file.close();
    t->nodes=vector<Node *>();
    if (line.back()==';') line.pop_back();
    auto root = new Node(line, &t->nodes);
    root->set_father(nullptr);
    t->n=(t->nodes.size()+1)/2;
    return t;
}

string Tree::to_string() {
    for (size_t i=n;i<nodes.size();i++) nodes[i]->set_name("NODE"+std::to_string(i));
    return nodes[n]->to_string() + ";";
}

Tree::~Tree() {
    for (size_t i = 0; i < nodes.size(); i++) {delete nodes[i];}
}

void Tree::write(const string &filename) {
    std::ofstream file;
    file.open(filename);
    if (file.is_open() == false) {
        cerr << "Error opening file " << filename << endl;
        return;
    }
    file << to_string() << endl;
    file.close();
}

bool Tree::check_structure() {
    if (n * 2 - 1 - nodes.size() != 0) return false;
    if (nodes[n]->get_father()!=nullptr) return false;
    for (size_t i = 0; i < nodes.size(); i++) {
        if (nodes[i]->get_left() != nullptr && nodes[i]->get_left()->get_father() != nodes[i]) return false;
        if (nodes[i]->get_right() != nullptr && nodes[i]->get_right()->get_father() != nodes[i]) return false;
    }
    return true;
}

int Tree::simulate_sequences(size_t l) {
    delete nodes[n]->get_seq();
    nodes[n]->set_seq(new vector<uint8_t>(l));
    std::uniform_int_distribution unif4(0, 3);
    for (size_t i = 0; i < l; i++) {
        nodes[n]->get_seq()->at(i) = Alignment::nuc_from_int(unif4(rng));
    }
    int m=0;
    m+=nodes[n]->get_left ()->simulate_sequences(nodes[n]->get_seq());
    m+=nodes[n]->get_right()->simulate_sequences(nodes[n]->get_seq());
    Alignment* a = Alignment::make_from_tree(this,true,true);//to perform pattern analysis
    set_sequences(a,true);
    return(m);
}

Tree* Tree::random_tree(size_t n,double avg_bra_len) {
    auto t = new Tree();
    t->n=n;
    t->nodes=vector<Node *>(n+n-1);
    for (size_t i = 0; i < 2*n-1; i++) {
        t->nodes[i]=new Node();
        if (i<n) t->nodes[i]->set_name("seq" + std::to_string(i + 1));
    }
    auto to_merge = vector<Node *>(n);
    std::exponential_distribution rexp(1.0/avg_bra_len);
    for (size_t i = 0; i < n; i++) {
        to_merge[i]=t->nodes[i];
        to_merge[i]->set_dist_father(rexp(rng));
    }
    size_t m=n+1;
    while (to_merge.size() > 2) {
        std::uniform_int_distribution unif(0, static_cast<int>(to_merge.size() - 1));
        int a = unif(rng);
        int b = a;
        while (a == b) b = unif(rng);
        t->nodes[m]->set_left(to_merge[a]);
        t->nodes[m]->set_right(to_merge[b]);
        to_merge[a]->set_father(t->nodes[m]);
        to_merge[b]->set_father(t->nodes[m]);
        t->nodes[m]->set_dist_father(rexp(rng));
        to_merge[a] = t->nodes[m];
        to_merge[b] = to_merge.back();
        to_merge.pop_back();
        m++;
    }
    t->nodes[n]->set_left (to_merge[0]);
    t->nodes[n]->set_right(to_merge[1]);
    to_merge[0]->set_father(t->nodes[n]);
    to_merge[1]->set_father(t->nodes[n]);
    return t;
}

Tree* Tree::nj(Alignment * ali) {
    size_t n=ali->get_n();
    auto t = new Tree();
    t->n=n;
    t->nodes=vector<Node *>(n+n-1);
    for (size_t i = 0; i < 2*n-1; i++) {
        t->nodes[i]=new Node();
        if (i<n) t->nodes[i]->set_name(ali->get_name(i));
    }
    auto to_merge = vector<int>(n);
    for (size_t i = 0; i < n; i++) to_merge[i]=i;
    vector<vector<double> > dist(2 * n, vector<double>(2 * n, 0.0));
    for (size_t i = 0; i < n; i++) for (size_t j = 0; j < i; j++) {
        dist[i][j] = ali->dist(i, j);
        dist[j][i] = dist[i][j];
    }
    vector<vector<double> > q(2 * n, vector<double>(2 * n, 0.0));
    vector<double> r(2 * n, 0.0);

    int m = 0;
    while (true) {
        for (size_t i = 0; i < to_merge.size(); i++) {
            r[to_merge[i]] = 0;
            for (size_t j = 0; j < to_merge.size(); j++) r[to_merge[i]] += dist[to_merge[i]][to_merge[j]];
        }
        if (to_merge.size() == 3) break;
        for (size_t i = 0; i < to_merge.size(); i++)
            for (size_t j = 0; j < i; j++)
                q[to_merge[i]][to_merge[j]] = (n - m - 2) * dist[to_merge[i]][to_merge[j]] - r[to_merge[i]] - r[to_merge[j]];
        size_t f = 1;
        size_t g = 0;
        for (size_t i = 0; i < to_merge.size(); i++)
            for (size_t j = 0; j < i; j++)
                if (q[to_merge[i]][to_merge[j]] < q[to_merge[f]][to_merge[g]]) {
                    f = i;
                    g = j;
                }
        Node *node = t->nodes[n+m];
        node->set_left(t->nodes[to_merge[f]]);
        node->set_right(t->nodes[to_merge[g]]);
        t->nodes[to_merge[f]]->set_father(node);
        t->nodes[to_merge[g]]->set_father(node);
        t->nodes[to_merge[f]]->set_dist_father(
            0.5 * dist[to_merge[f]][to_merge[g]] + 0.5 / (n - m - 2) * (r[to_merge[f]] - r[to_merge[g]]));
        t->nodes[to_merge[g]]->set_dist_father(dist[to_merge[f]][to_merge[g]] - t->nodes[to_merge[f]]->get_dist_father());
        for (size_t k = 0; k < to_merge.size(); k++) {
            dist[n + m][to_merge[k]] = 0.5 * (dist[to_merge[f]][to_merge[k]] + dist[to_merge[g]][to_merge[k]] - dist[
                                                  to_merge[f]][to_merge[g]]);
            dist[to_merge[k]][n + m] = dist[n + m][to_merge[k]];
        }
        to_merge[f] = n + m;
        to_merge[g] = to_merge.back();
        to_merge.pop_back();
        m++;
    }
    t->nodes[to_merge[0]]->set_dist_father(
        0.5 * dist[to_merge[0]][to_merge[1]] + 0.5 / (n - m - 2) * (r[to_merge[0]] - r[to_merge[1]]));
    t->nodes[to_merge[1]]->set_dist_father(dist[to_merge[0]][to_merge[1]] - t->nodes[to_merge[0]]->get_dist_father());
    t->nodes[to_merge[2]]->set_dist_father(
        0.5 * (dist[to_merge[0]][to_merge[2]] + dist[to_merge[1]][to_merge[2]] - dist[to_merge[0]][to_merge[1]]));
    t->nodes[n+m]->set_left (t->nodes[to_merge[0]]);
    t->nodes[n+m]->set_right(t->nodes[to_merge[1]]);
    t->nodes[to_merge[0]]->set_father(t->nodes[n+m]);
    t->nodes[to_merge[1]]->set_father(t->nodes[n+m]);
    t->nodes[n+m+1]->set_left (t->nodes[to_merge[2]]);
    t->nodes[n+m+1]->set_right(t->nodes[n+m]);
    t->nodes[to_merge[2]]->set_father(t->nodes[n+m+1]);
    t->nodes[n+m]->set_father(t->nodes[n+m+1]);
    t->nodes[n+m]->set_dist_father(0.0);
    t->nodes[n+m]->set_dist_father(t->nodes[to_merge[2]]->get_dist_father()*0.5);
    t->nodes[to_merge[2]]->set_dist_father(t->nodes[to_merge[2]]->get_dist_father()*0.5);
    std::swap(t->nodes[n+m+1],t->nodes[n]);
    t->nodes[n]->set_father(nullptr);
    //for (i=0;i<nodes.size();i++) if (nodes[i]->get_dist_father()<=0.0) nodes[i]->set_dist_father(1e-5);
    return t;
}

void Tree::set_sequences(Alignment *ali,bool is_ordered) {
    for (size_t i=n;i<nodes.size();i++) {
        delete nodes[i]->get_seq();
        nodes[i]->set_seq(new vector<uint8_t>(ali->get_npat(),0b00010000));
    }
    if (is_ordered) {
        for (size_t i=0;i<ali->get_n();i++) get_node(i)->set_seq(ali->get_seq(i));
    } else {
        for (size_t i=0;i<ali->get_n();i++) for (size_t j=0;j<nodes.size();j++)
            if (get_node(j)->get_name()==ali->get_name(i)) get_node(j)->set_seq(ali->get_seq(i));
    }
    npat=ali->get_npat();
    patindex=ali->get_patindex();
    weights=ali->get_weights();
}

int Tree::parsimony() {
    //Determine order
    vector<Node*> todo;
    todo.push_back(nodes[n]);
    size_t i=0;
    while (i<todo.size()) {
        if (todo[i]->get_left()!=nullptr) {
            todo.push_back(todo[i]->get_left());
            todo.push_back(todo[i]->get_right());
        }
        i++;
    }

    //Fitch first pass
    int score=0;
    for (i=todo.size()-1;true;i--) {
        for (size_t pos=0;pos<npat;pos++) {
            if (todo[i]->get_left()==nullptr) {
                if (Alignment::nuc_is_unsure(todo[i]->get_seq()->at(pos))) todo[i]->get_seq()->at(pos)=0b00011111;
            } else {
                uint8_t toset=todo[i]->get_left()->get_seq()->at(pos)&todo[i]->get_right()->get_seq()->at(pos)&0b00001111;//intersect
                if (toset==0u) {
                    score+=weights->at(pos);
                    toset=todo[i]->get_left()->get_seq()->at(pos)|todo[i]->get_right()->get_seq()->at(pos);//union
                }
                Alignment::nuc_set_as(&(todo[i]->get_seq()->at(pos)),toset);//preserve missing status
            }
        }
        if (i==0) break;
    }

    //Fitch second pass
    for (i=0;i<todo.size();i++) for (size_t pos=0;pos<npat;pos++) {
        uint8_t toset;
        if (todo[i]->get_father()==nullptr) toset=todo[i]->get_seq()->at(pos)&0b00001111;
        else {
            toset=todo[i]->get_seq()->at(pos)&todo[i]->get_father()->get_seq()->at(pos)&0b00001111;
            if (toset==0u) toset=todo[i]->get_seq()->at(pos)&0b00001111;
        }
        toset=toset&-toset;//pick first common bit
        Alignment::nuc_set_as(&(todo[i]->get_seq()->at(pos)),toset);//preserve missing status

    }
    return(score);
}

int Tree::parsimony_recursive() {
    int score=0;
    for (size_t pos=0;pos<npat;pos++) {
        score+=nodes[n]->fitch_pass1(pos)*weights->at(pos);
        nodes[n]->fitch_pass2(pos);
    }
    return(score);
}

int Tree::count_subs() {
    int subs=0;
    for (size_t pos=0;pos<npat;pos++) {
        for (size_t i=0;i<nodes.size();i++) {
            if (nodes[i]->get_father()==nullptr) continue;
            if (!Alignment::nuc_is_currently_same(nodes[i]->get_seq()->at(pos),nodes[i]->get_father()->get_seq()->at(pos))) subs+=weights->at(pos);
        }
    }
    return subs;
}

double Tree::sum_bra_len() {
    double res=0;
        for (size_t i=0;i<nodes.size();i++) {
            if (nodes[i]->get_father()==nullptr) continue;
            res+=nodes[i]->get_dist_father();
        }
    return res;
}

double Tree::likelihood_anc() {
    double ll=0.0;
    double len,ex,pid,pnid;
    int s1,s2;
    for (size_t i=0;i<nodes.size();i++) {
        if (nodes[i]->get_father()==nullptr) continue;
        len=nodes[i]->get_dist_father();
        ex=exp(-len);
        pid =log(0.25+0.75*ex);
        pnid=log(0.25-0.25*ex);
        s1=0;
        s2=0;
            for (size_t j=0;j<npat;j++) {
                auto seq1=nodes[i]->get_father()->get_seq();
                auto seq2=nodes[i]->get_seq();
                if (Alignment::nuc_is_currently_same(seq1->at(j),seq2->at(j)))
                    s1+=weights->at(j);
                else s2+=weights->at(j);
        }
        ll+=pid*s1+pnid*s2;
    }
    ll=ll+log(0.25)*patindex->size();
    return ll;
}
