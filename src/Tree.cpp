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
    std::swap(t->nodes[n*2-2],t->nodes[n]);
    t->reorder();
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
    std::swap(t->nodes[t->n*2-2],t->nodes[t->n]);
    t->reorder();
    return t;
}

string Tree::to_string() {
    for (size_t i=n;i<nodes.size();i++) nodes[i]->set_name("NODE"+std::to_string(i));
    return nodes[2*n-2]->to_string() + ";";
}

Tree::~Tree() {
    for (size_t i = n; i < nodes.size(); i++) {delete nodes[i]->get_seq();}
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
    if (nodes[2*n-2]->get_father()!=nullptr) return false;
    for (size_t i = 0; i < nodes.size(); i++) {
        if (nodes[i]->get_left() != nullptr && nodes[i]->get_left()->get_father() != nodes[i]) return false;
        if (nodes[i]->get_right() != nullptr && nodes[i]->get_right()->get_father() != nodes[i]) return false;
        if (nodes[i]->get_index() != i) return false;
    }
    return true;
}

int Tree::simulate_sequences(size_t l) {
    for (int i=0;i<2*(int)n-1;i++) delete nodes[i]->get_seq();
    nodes[2*n-2]->set_seq(new vector<uint8_t>(l));
    std::uniform_int_distribution unif4(0, 3);
    std::uniform_int_distribution unif_l(0,static_cast<int>(l)-1);
    for (size_t i = 0; i < l; i++) {
        nodes[2*n-2]->get_seq()->at(i) = Alignment::nuc_from_int(unif4(rng));
    }
    int m=0;
    for (int z=2*(int)n-3;z>=0;z--) {
        vector<uint8_t>*above=nodes[z]->get_father()->get_seq();
        nodes[z]->set_seq(new vector<uint8_t>(*above));
        std::poisson_distribution pois(nodes[z]->get_dist_father()*l);
        int m2=pois(rng);
        for (int i=0;i<m2;i++) {
            int p=unif_l(rng);
            while (above->at(p)==nodes[z]->get_seq()->at(p)) {
                nodes[z]->get_seq()->at(p)=Alignment::nuc_from_int(unif4(rng));
            }
        }
        m+=m2;
    }
    Alignment* a = Alignment::make_from_tree(this,true,true);//to perform pattern analysis
    for (size_t i=0;i<nodes.size();i++) nodes[i]->set_seq(nullptr);
    set_sequences(a);
    a->unlink();
    delete(a);
    return(m);
}

int Tree::simulate_sequences_recursive(size_t l) {
    delete nodes[2*n-2]->get_seq();
    nodes[2*n-2]->set_seq(new vector<uint8_t>(l));
    std::uniform_int_distribution unif4(0, 3);
    for (size_t i = 0; i < l; i++) {
        nodes[2*n-2]->get_seq()->at(i) = Alignment::nuc_from_int(unif4(rng));
    }
    int m=0;
    m+=nodes[2*n-2]->get_left ()->simulate_sequences(nodes[2*n-2]->get_seq());
    m+=nodes[2*n-2]->get_right()->simulate_sequences(nodes[2*n-2]->get_seq());
    Alignment* a = Alignment::make_from_tree(this,true,true);//to perform pattern analysis
    for (size_t i=0;i<nodes.size();i++) nodes[i]->set_seq(nullptr);
    set_sequences(a);
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
    std::swap(t->nodes[n*2-2],t->nodes[n]);
    t->reorder();
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
    auto to_merge = vector<size_t>(n);
    for (size_t i = 0; i < n; i++) to_merge[i]=i;
    vector<vector<double> > dist(2 * n, vector<double>(2 * n, 0.0));
    for (size_t i = 0; i < n; i++) for (size_t j = 0; j < i; j++) {
        dist[i][j] = ali->dist(i, j);
        dist[j][i] = dist[i][j];
    }
    vector<vector<double> > q(2 * n, vector<double>(2 * n, 0.0));
    vector<double> r(2 * n, 0.0);

    size_t m = 0;
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
    t->nodes[n+m+1]->set_father(nullptr);
    //for (i=0;i<nodes.size();i++) if (nodes[i]->get_dist_father()<=0.0) nodes[i]->set_dist_father(1e-5);
    t->reorder();
    return t;
}

void Tree::set_sequences(Alignment *ali) {
    for (size_t i=0;i<nodes.size();i++) {
        delete nodes[i]->get_seq();
    }
    for (size_t i=0;i<ali->get_n();i++) for (size_t j=0;j<nodes.size();j++)
        if (get_node(j)->get_name()==ali->get_name(i)) get_node(j)->set_seq(ali->get_seq(i));
    for (size_t i=0;i<nodes.size();i++) {
        if (nodes[i]->get_seq()==nullptr) nodes[i]->set_seq(new vector<uint8_t>(ali->get_npat(),0b00010000));
    }
    npat=ali->get_npat();
    patindex=ali->get_patindex();
    weights=ali->get_weights();
}

int Tree::parsimony() {
    //Fitch first pass
    int i=0;
    int score=0;
    for (i=0;i<(int)nodes.size();i++) {
        for (size_t pos=0;pos<npat;pos++) {
            if (nodes[i]->get_left()==nullptr) {
                if (Alignment::nuc_is_unsure(nodes[i]->get_seq()->at(pos))) nodes[i]->get_seq()->at(pos)=0b00011111;
            } else {
                uint8_t toset=nodes[i]->get_left()->get_seq()->at(pos)&nodes[i]->get_right()->get_seq()->at(pos)&0b00001111;//intersect
                if (toset==0u) {
                    score+=weights->at(pos);
                    toset=nodes[i]->get_left()->get_seq()->at(pos)|nodes[i]->get_right()->get_seq()->at(pos);//union
                }
                Alignment::nuc_set_as(&(nodes[i]->get_seq()->at(pos)),toset);//preserve missing status
            }
        }
    }

    //Fitch second pass
    for (i=nodes.size()-1;i>=0;i--) {
        for (size_t pos=0;pos<npat;pos++) {
            uint8_t toset;
            if (nodes[i]->get_father()==nullptr) toset=nodes[i]->get_seq()->at(pos)&0b00001111;
            else {
                toset=nodes[i]->get_seq()->at(pos)&nodes[i]->get_father()->get_seq()->at(pos)&0b00001111;
                if (toset==0u) toset=nodes[i]->get_seq()->at(pos)&0b00001111;
            }
            toset=toset&-toset;//pick first common bit
            Alignment::nuc_set_as(&(nodes[i]->get_seq()->at(pos)),toset);//preserve missing status
        }
    }
    return(score);
}

int Tree::parsimony_recursive() {
    int score=0;
    for (size_t pos=0;pos<npat;pos++) {
        score+=nodes[2*n-2]->fitch_pass1(pos)*weights->at(pos);
        nodes[2*n-2]->fitch_pass2(pos);
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
        ex = exp(-4.0*len/3.0);
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

void Tree::reorder() {
    vector<Node*> todo;
    todo.push_back(nodes[2*n-2]);
    size_t i=0,j=0;
    while (i<todo.size()) {
        if (todo[i]->get_left()!=nullptr) {
            todo.push_back(todo[i]->get_left());
            todo.push_back(todo[i]->get_right());
        }
        i++;
    }
    for (i=0;i<2*n-1;i++) {
        if (todo[2*n-2-i]->get_left()!=nullptr) {
            nodes[n+j]=todo[2*n-2-i];
            j++;
        }
    }
    for (i=0;i<2*n-1;i++) nodes[i]->set_index(i);
}

double Tree::jointml() {
    vector<vector<double> > L(nodes.size(), vector<double>(4, 0.0));
    vector<vector<uint8_t> > C(nodes.size(), vector<uint8_t>(4, 0));
    vector<double> tmp(4, 0.0);
    double res = 0.0;
    int z,w;

    vector<double> pid(nodes.size(),0);
    vector<double> pnid(nodes.size(),0);
    for (z = 0; z < (int)n*2-1; z++) {
        size_t zind=nodes[z]->get_index();
        double ex = exp(-4.0*nodes[z]->get_dist_father()/3.0);
        pid[zind] = log(0.25 + 0.75 * ex);
        pnid[zind] = log(0.25 - 0.25 * ex);
    }

    for (size_t pos = 0; pos < npat; pos++) {
        //First pass from leaves to root (excluded)
        for (z = 0; z < (int)n*2-1; z++) {
            size_t zind=nodes[z]->get_index();

            for (int i = 0; i < 4; i++) {
                if (z<(int)n) {
                    if (Alignment::nuc_is_unsure(nodes[z]->get_seq()->at(pos))) {
                        L[zind][i] = pid[zind];
                        C[zind][i] = 1<<i;
                    } else {
                        C[zind][i] = nodes[z]->get_seq()->at(pos);
                        if (i == __builtin_ctz(nodes[z]->get_seq()->at(pos)))
                            L[zind][i] = pid[zind];
                        else
                            L[zind][i] = pnid[zind];
                    }
                } else {
                    size_t left=nodes[z]->get_left()->get_index();
                    size_t right=nodes[z]->get_right()->get_index();
                    for (int j=0;j<4;j++) {
                        tmp[j]=L[left][j]+L[right][j];
                        if (i==j) tmp[j]+=pid[zind]; else tmp[j]+=pnid[zind];
                    }
                    w=0;
                    for (size_t j=1;j<4;j++) {if (tmp[j]>tmp[w]) w=j;}
                    L[zind][i]=tmp[w];
                    C[zind][i]=1<<w;
                }
            }
        }

        //Root
        z=2*n-2;
        size_t left=nodes[z]->get_left()->get_index();
        size_t right=nodes[z]->get_right()->get_index();
        for (size_t j=0;j<4;j++) {
            L[z][j]=L[left][j]+L[right][j]+log(0.25);
        }
        w=0;
        for (size_t j=1;j<4;j++) {if (L[z][j]>L[z][w]) w=j;}
        res+=L[z][w]*weights->at(pos);
        Alignment::nuc_set_as(&(nodes[z]->get_seq()->at(pos)), 1<<w); //preserve missing status

        //Second pass from root (excluded) to leaves
        for (z = 2*n-3;z>=0; z--) {
            size_t zind=nodes[z]->get_index();
            size_t i=__builtin_ctz(nodes[z]->get_father()->get_seq()->at(pos)&0b00001111);
            Alignment::nuc_set_as(&(nodes[z]->get_seq()->at(pos)), C[zind][i]); //preserve missing status
        }
    }
    return (res);
}

double Tree::likelihood() {
    vector<vector<double> > L(nodes.size(), vector<double>(4, 0.0));
    double res = 0.0;

    vector<double> pid(nodes.size(),0);
    vector<double> pnid(nodes.size(),0);
    for (size_t z = 0; z < nodes.size()-1; z++) {
        size_t zind=nodes[z]->get_index();
        double ex = exp(-4.0*nodes[z]->get_dist_father()/3.0);
        pid[zind] = (0.25 + 0.75 * ex);
        pnid[zind] = (0.25 - 0.25 * ex);
    }

    for (size_t pos = 0; pos < npat; pos++) {
        //Pass from leaves to root (included)
        for (size_t z = 0; z < nodes.size(); z++) {
            size_t zind=nodes[z]->get_index();
            for (int a = 0; a < 4; a++) {
                if (nodes[z]->get_left() == nullptr) {
                    if (Alignment::nuc_is_unsure(nodes[z]->get_seq()->at(pos))) {
                        L[zind][a] = 1;
                        cout<<"not happening"<<endl;
                    } else {
                        if (a == __builtin_ctz(nodes[z]->get_seq()->at(pos)))
                            L[zind][a] = 1;
                        else
                            L[zind][a] = 0;
                    }
                } else {
                    size_t left=nodes[z]->get_left()->get_index();
                    size_t right=nodes[z]->get_right()->get_index();
                    double l=0,r=0;
                    for (int b=0;b<4;b++) {
                        l+=L[left][b]*(a==b?pid[left]:pnid[left]);
                    }
                    for (int c=0;c<4;c++) {
                        r+=L[right][c]*(a==c?pid[right]:pnid[right]);
                    }
                    L[zind][a]=l*r;
                }
            }
        }
        size_t root=nodes.size()-1;
        res+=log(0.25*(L[root][0]+L[root][1]+L[root][2]+L[root][3]))*weights->at(pos);
    }
    return res;
}