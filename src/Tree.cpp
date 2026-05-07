//
// Created by Didelot, Xavier on 29/09/2025.

#include "Tree.h"
#include "Alignment.h"

#include <iostream>

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
    if (verbose) vprint("Opening file {}",filename);
    file.open(filename);
    if (file.is_open() == false) {
        std::println(std::cerr,"Error opening file {}",filename);
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
    Node * root=nodes[2*n-2];
    return "("+root->get_left()->get_left()->to_string() +","+root->get_left()->get_right()->to_string()+","+root->get_right()->to_string() + ")"+root->get_name()+";";
}

Tree::~Tree() {
    for (size_t i = 0; i < nodes.size(); i++) {delete nodes[i];}
}

void Tree::write(const string &filename) {
    std::ofstream file;
    file.open(filename);
    if (file.is_open() == false) {
        std::println(std::cerr,"Error opening file {}",filename);
        return;
    }
    std::println(file,"{}",to_string());
    file.close();
}

bool Tree::check_structure() {
    if (n * 2 - 1 - nodes.size() != 0) {
        std::println(std::cerr,"Incorrect number of nodes");
        return false;
    }
    if (nodes[2*n-2]->get_father()!=nullptr) {
        std::println(std::cerr,"Root is not last");
        return false;
    }
    if (nodes[2*n-2]->get_left()->get_left()==nullptr) {
        std::println(std::cerr,"Invalid root");
        return false;
    }
    if (nodes[2*n-2]->get_left()->get_dist_father()!=0.0) {
        std::println(std::cerr,"Invalid root distance");
        return false;
    }
    if (nodes[2*n-3]!=nodes[2*n-2]->get_left()) {
        std::println(std::cerr,"Null node is not penultimate");
        return false;
    }
    for (size_t i = 0; i < nodes.size(); i++) {
        if (i<nodes.size()-2 && nodes[i]->get_dist_father()==0.0) {
            std::println(std::cerr,"Node with distance zero other than root or null node");
            return false;
        }
        if (nodes[i]->get_left() != nullptr && nodes[i]->get_left()->get_father() != nodes[i]) {
            std::println(std::cerr,"Incorrect left pointer");
            return false;
        }
        if (nodes[i]->get_right() != nullptr && nodes[i]->get_right()->get_father() != nodes[i]) {
            std::println(std::cerr,"Incorrect right pointer");
            return false;
        }
        if (nodes[i]->get_index() != i) {
            std::println(std::cerr,"Incorrect node index");
            return false;
        }
    }
    return true;
}

int Tree::simulate_sequences(size_t l) {
    for (int i=0;i<2*(int)n-1;i++) delete nodes[i]->get_seq();
    nodes[2*n-2]->set_seq(new vector<uint8_t>(l));
    std::uniform_int_distribution unif_l(0,static_cast<int>(l)-1);
    for (size_t i = 0; i < l; i++) {
        nodes[2*n-2]->get_seq()->at(i) = Alignment::random_nuc();
    }
    int m=0;
    for (int z=2*(int)n-3;z>=0;z--) {
        vector<uint8_t>*above=nodes[z]->get_father()->get_seq();
        nodes[z]->set_seq(new vector<uint8_t>(*above));
        int m2=0;
        if (nodes[z]->get_dist_father()>0.0) {
            std::poisson_distribution pois(nodes[z]->get_dist_father()*l);
            m2=pois(rng);
        }
        for (int i=0;i<m2;i++) {
            int p=unif_l(rng);
            nodes[z]->get_seq()->at(p)=Alignment::random_nuc_diff(above->at(p));
        }
        m+=m2;
    }
    simulate_pattern_analysis();
    return(m);
}

int Tree::simulate_sequences_rec(size_t l,double rho,double delta,double nu) {
    for (int i=0;i<2*(int)n-1;i++) delete nodes[i]->get_seq();
    nodes[2*n-2]->set_seq(new vector<uint8_t>(l));
    std::uniform_real_distribution runif(0.0,1.0);
    for (size_t i = 0; i < l; i++) {
        nodes[2*n-2]->get_seq()->at(i) = Alignment::random_nuc();
    }
    int m=0;
    for (int z=2*(int)n-3;z>=0;z--) {
        vector<uint8_t>*above=nodes[z]->get_father()->get_seq();
        nodes[z]->set_seq(new vector<uint8_t>(*above));
        double a=1.0-exp(-rho*nodes[z]->get_dist_father());
        double b=1.0-exp(-1.0/delta);
        bool rec=runif(rng)<a/(a+b);
        size_t curpos=0;
        while (curpos<l) {
            size_t upto;
            int m2;
            if (nodes[z]->get_dist_father()>0) {
                std::geometric_distribution geo(rec?b:a);
                upto=curpos+geo(rng)+1;
                if (upto>l) upto=l;
                std::poisson_distribution pois((rec?nu:nodes[z]->get_dist_father())*(upto-curpos));
                m2=pois(rng);
            } else {
                upto=l;
                m2=0;
            }
            for (int i=0;i<m2;i++) {
                int p=curpos+(rng()%(upto-curpos));
                nodes[z]->get_seq()->at(p)=Alignment::random_nuc_diff(above->at(p));
            }
            m+=m2;
            rec=!rec;
            curpos=upto;
        }
    }
    simulate_pattern_analysis();
    return(m);
}

void Tree::simulate_pattern_analysis() {
    size_t l=nodes[0]->get_seq()->size();
    Alignment* a  = Alignment::make_from_tree(this,true);
    Alignment* a2 = Alignment::make_from_tree(this,false);a2->pattern_analysis();
    vector<size_t> anti(a2->get_npat());
    for (int j=static_cast<int>(l)-1;j>=0;j--) anti[a2->get_patindex()->at(j)]=j;
    for (size_t i=0;i<nodes.size();i++) {
        delete nodes[i]->get_seq();
        nodes[i]->set_seq(new vector<uint8_t>(a2->get_npat()));
        for (size_t j=0;j<a2->get_npat();j++) nodes[i]->get_seq()->at(j)=a->get_data(i,anti[j]);
    }
    delete(a);
    ali=a2;
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

Tree* Tree::nj(vector<vector<double>> * d,vector<string> * names) {
    size_t n=names->size();
    auto t = new Tree();
    t->n=n;
    t->nodes=vector<Node *>(n+n-1);
    for (size_t i = 0; i < 2*n-1; i++) {
        t->nodes[i]=new Node();
        if (i<n) t->nodes[i]->set_name(names->at(i));
    }
    auto to_merge = vector<size_t>(n);
    for (size_t i = 0; i < n; i++) to_merge[i]=i;
    vector<vector<double> > dist(2 * n, vector<double>(2 * n, 0.0));
    for (size_t i = 0; i < n; i++) for (size_t j = 0; j < i; j++) {
        dist[i][j] = (*d)[i][j];
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
    t->nodes[n+m]->set_dist_father(t->nodes[to_merge[2]]->get_dist_father());
    t->nodes[to_merge[2]]->set_dist_father(0.0);
    t->nodes[n+m+1]->set_father(nullptr);
    t->reorder();
    for (size_t i=0;i<t->nodes.size()-1;i++)
        if (t->nodes[i]!=t->nodes[2*n-2]->get_left() && t->nodes[i]->get_dist_father()<=1e-6) t->nodes[i]->set_dist_father(1e-6);
    return t;
}

void Tree::set_sequences(Alignment *ali) {
    for (size_t i=0;i<nodes.size();i++) {
        delete nodes[i]->get_seq();
    }
    for (size_t i=0;i<ali->get_n();i++) for (size_t j=0;j<nodes.size();j++)
        if (get_node(j)->get_name()==ali->get_name(i)) get_node(j)->set_seq(new vector(*ali->get_seq(i)));
    for (size_t i=0;i<nodes.size();i++) {
        if (nodes[i]->get_seq()==nullptr) nodes[i]->set_seq(new vector<uint8_t>(ali->get_npat(),0b00010000));
    }
    this->ali=ali;
}

int Tree::parsimony() {
    //Fitch first pass
    int i=0;
    int score=0;
    for (i=0;i<(int)nodes.size();i++) {
        for (size_t pos=0;pos<ali->get_npat();pos++) {
            if (nodes[i]->get_left()==nullptr) {
                if (Alignment::nuc_is_unsure(nodes[i]->get_seq()->at(pos))) nodes[i]->get_seq()->at(pos)=0b00011111;
            } else {
                uint8_t toset=nodes[i]->get_left()->get_seq()->at(pos)&nodes[i]->get_right()->get_seq()->at(pos)&0b00001111;//intersect
                if (toset==0u) {
                    score+=ali->get_weights()->at(pos);
                    toset=nodes[i]->get_left()->get_seq()->at(pos)|nodes[i]->get_right()->get_seq()->at(pos);//union
                    if (nodes[i]->get_father()==nullptr) toset=nodes[i]->get_left()->get_seq()->at(pos);//special case for root since working with unrooted trees
                }
                Alignment::nuc_set_as(&(nodes[i]->get_seq()->at(pos)),toset);//preserve missing status
            }
        }
    }

    //Fitch second pass
    for (i=nodes.size()-1;i>=0;i--) {
        for (size_t pos=0;pos<ali->get_npat();pos++) {
            uint8_t toset;
            if (nodes[i]->get_father()==nullptr) toset=nodes[i]->get_seq()->at(pos)&0b00001111;
            else {
                toset=nodes[i]->get_seq()->at(pos)&nodes[i]->get_father()->get_seq()->at(pos)&0b00001111;
                if (toset==0u && nodes[i]->get_dist_father()==0) std::println("Problem: adding a substitution on a branch of length zero!");
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
    for (size_t pos=0;pos<ali->get_npat();pos++) {
        score+=nodes[2*n-2]->fitch_pass1(pos)*ali->get_weights()->at(pos);
        nodes[2*n-2]->fitch_pass2(pos);
    }
    return(score);
}

int Tree::count_subs() {
    int subs=0;
    for (size_t pos=0;pos<ali->get_npat();pos++) {
        for (size_t i=0;i<nodes.size();i++) {
            if (nodes[i]->get_father()==nullptr) continue;
            if (!Alignment::nuc_is_currently_same(nodes[i]->get_seq()->at(pos),nodes[i]->get_father()->get_seq()->at(pos))) subs+=ali->get_weights()->at(pos);
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
    vector<uint8_t> * seq1;
    vector<uint8_t> * seq2;
    for (size_t i=0;i<nodes.size();i++) {
        if (nodes[i]->get_father()==nullptr) continue;
        seq1=nodes[i]->get_father()->get_seq();
        seq2=nodes[i]->get_seq();
        len=nodes[i]->get_dist_father();
        ex = exp(-4.0*len/3.0);
        pid =log(0.25+0.75*ex);
        pnid=log(0.25-0.25*ex);
        s1=0;
        s2=0;
            for (size_t j=0;j<ali->get_npat();j++) {
                if (Alignment::nuc_is_currently_same(seq1->at(j),seq2->at(j)))
                    s1+=ali->get_weights()->at(j);
                else s2+=ali->get_weights()->at(j);
        }
        if (s2>0  && len==0.0) std::println("Problem: there are substitutions on a branch of length zero!");
        if (s2==0 && len==0.0) pnid=0.0;//otherwise get nan below
        ll+=pid*s1+pnid*s2;
    }
    ll=ll+log(0.25)*ali->get_patindex()->size();
    return ll;
}

double Tree::likelihood_anc_rec(double rho, double delta, double nu) {
    if (profiling) Profiler::prof.start(__func__);
    double ll=0.0;
    size_t L=ali->get_patindex()->size();
    #pragma omp parallel for reduction(+:ll) schedule(static)
    for (size_t i=0;i<nodes.size();i++) {
        if (nodes[i]->get_father()==nullptr) continue;
        vector<uint8_t> * seq1=nodes[i]->get_father()->get_seq();
        vector<uint8_t> * seq2=nodes[i]->get_seq();
        double len=nodes[i]->get_dist_father();
        double ex = exp(-4.0*len/3.0);
        double pid =0.25+0.75*ex;
        double pnid=0.25-0.25*ex;
        double a=1.0-exp(-rho*len);
        //double b=1.0-exp(-1.0/delta);
        double b;
        if (a==0.0) b=0.5; else b=1.0-exp(-rho*len/expm1(rho*len/delta));//Using M/M/inf queue
        double alpha1=b/(a+b);
        double alpha2=a/(a+b);
        vector gamma={1.0-a, a, b, 1.0-b};
        size_t maxskip=100;
        vector gammas(maxskip, vector(4, 0.0));
        gammas[0]=gamma;
        for (size_t k=1;k<maxskip;k++) {
            gammas[k][0]=gammas[k-1][0]*pid*gamma[0]+gammas[k-1][1]*(1.0-nu)*gamma[2];
            gammas[k][1]=gammas[k-1][0]*pid*gamma[1]+gammas[k-1][1]*(1.0-nu)*gamma[3];
            gammas[k][2]=gammas[k-1][2]*pid*gamma[0]+gammas[k-1][3]*(1.0-nu)*gamma[2];
            gammas[k][3]=gammas[k-1][2]*pid*gamma[1]+gammas[k-1][3]*(1.0-nu)*gamma[3];
        }
        double lscale=0.0;
        size_t j=0;
        while (j<L) {//Forward algorithm
            size_t skipping=0;
            bool is_same=true;
            while (is_same) {
                size_t pat=(*ali->get_patindex())[j+skipping];
                is_same=Alignment::nuc_is_currently_same((*seq1)[pat],(*seq2)[pat]);
                if (j+skipping>L-2 || skipping>maxskip-2) break;
                if (is_same) skipping++;
            }
            double p1,p2;
            if (is_same) {p1=pid;p2=1.0-nu;} else {p1=pnid;p2=nu;}
            double newalpha1=(alpha1*gammas[skipping][0]+alpha2*gammas[skipping][2])*p1;
            alpha2=          (alpha1*gammas[skipping][1]+alpha2*gammas[skipping][3])*p2;
            alpha1=newalpha1;
            double scale=alpha1+alpha2;
            alpha1=alpha1/scale;
            alpha2=alpha2/scale;
            lscale+=log(scale);
            j+=skipping+1;
        }
        ll+=lscale;
    }
    ll=ll+log(0.25)*static_cast<double>(L);
    if (profiling) Profiler::prof.stop(__func__);
    return ll;
}

void Tree::reorder() {
    if (profiling) Profiler::prof.start(__func__);
    vector<Node*> todo;
    //Find the root
    int a=0;
    while (nodes[2*n-2-a]->get_father()!=nullptr) a++;
    todo.push_back(nodes[2*n-2-a]);
    //If only a single leaf on the left of the root, swap left and right
    if (todo[0]->get_left()->get_left()==nullptr) {
        Node * n=todo[0]->get_left();
        todo[0]->set_left(todo[0]->get_right());
        todo[0]->set_right(n);
    }
    //Make the left branch of length zero
    double l=todo[0]->get_left()->get_dist_father();
    todo[0]->get_right()->set_dist_father(todo[0]->get_right()->get_dist_father()+l);
    todo[0]->get_left()->set_dist_father(0.0);
    //Order the nodes
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
    if (profiling) Profiler::prof.stop(__func__);

}

vector<vector<double>> * Tree::distmat() {
    if (profiling) Profiler::prof.start(__func__);
    auto d=new vector<vector<double>>(n, vector<double>(n, 0.0));
    for (size_t i = 0; i < n; i++) for (size_t j = 0; j < i; j++) {
        vector<Node*> path1;
        path1.push_back(nodes[i]);
        while (path1.back()->get_father()!=nullptr) path1.push_back(path1.back()->get_father());
        vector<Node*> path2;
        path2.push_back(nodes[j]);
        while (path2.back()->get_father()!=nullptr) {
            path2.push_back(path2.back()->get_father());
            if (std::find(path1.begin(), path1.end(), path2.back()) != path1.end()) break;
        }
        double di=0;
        Node * cur=nodes[i];
        while (cur!=path2.back()) {
            di+=cur->get_dist_father();
            cur=cur->get_father();
        }
        cur=nodes[j];
        while (cur!=path2.back()) {
            di+=cur->get_dist_father();
            cur=cur->get_father();
        }
        (*d)[i][j] = di;
        (*d)[j][i] = (*d)[i][j];
    }
    if (profiling) Profiler::prof.stop(__func__);
    return d;
}

double Tree::jointml() {
    if (profiling) Profiler::prof.start(__func__);
    vector<vector<double> > L(nodes.size(), vector<double>(4, 0.0));
    vector<vector<uint8_t> > C(nodes.size(), vector<uint8_t>(4, 0));
    double tmp[4];
    double res = 0.0;
    int z,w;
    const double log025 = log(0.25);

    vector<double> pid(nodes.size(),0);
    vector<double> pnid(nodes.size(),0);
    for (z = 0; z < (int)n*2-1; z++) {
        size_t zind=nodes[z]->get_index();
        double ex = exp(-4.0*nodes[z]->get_dist_father()/3.0);
        pid[zind] = log(0.25 + 0.75 * ex);
        pnid[zind] = log(0.25 - 0.25 * ex);
    }

    for (size_t pos = 0; pos < ali->get_npat(); pos++) {
        //First pass for leaves
        for (z = 0; z < (int)n; z++) {
            size_t zind=nodes[z]->get_index();
            uint8_t nuc=nodes[z]->get_seq()->at(pos);
            int nucint = __builtin_ctz(nuc);
            bool nuc_is_unsure =Alignment::nuc_is_unsure(nuc);
            for (int i = 0; i < 4; i++) {
                if (nuc_is_unsure) {
                    L[zind][i] = pid[zind];
                    C[zind][i] = 1<<i;
                } else {
                    C[zind][i] = nuc;
                    if (i == nucint)
                        L[zind][i] = pid[zind];
                    else
                        L[zind][i] = pnid[zind];
                }
            }
        }
        //First pass for internal nodes other than root
        for (z = (int)n; z < (int)n*2-1; z++) {
            size_t zind=nodes[z]->get_index();
            size_t left=nodes[z]->get_left()->get_index();
            size_t right=nodes[z]->get_right()->get_index();
            double base[4];
            for (int j = 0; j < 4; j++) base[j] = L[left][j] + L[right][j];
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    tmp[j] = base[j] + (i == j ? pid[zind] : pnid[zind]);
                }
                w = 0;
                for (int j = 1; j < 4; j++) if (tmp[j] > tmp[w]) w = j;
                L[zind][i] = tmp[w];
                C[zind][i] = 1 << w;
            }
        }

        //Root
        z=2*n-2;
        size_t left=nodes[z]->get_left()->get_index();
        size_t right=nodes[z]->get_right()->get_index();
        for (size_t j=0;j<4;j++) {
            L[z][j]=L[left][j]+L[right][j]+log025;
        }
        w=0;
        for (size_t j=1;j<4;j++) {if (L[z][j]>L[z][w]) w=j;}
        res+=L[z][w]*ali->get_weights()->at(pos);
        Alignment::nuc_set_as(&(nodes[z]->get_seq()->at(pos)), 1<<w); //preserve missing status

        //Second pass from root (excluded) to leaves
        for (z = 2*n-3;z>=0; z--) {
            size_t zind=nodes[z]->get_index();
            size_t i=__builtin_ctz(nodes[z]->get_father()->get_seq()->at(pos)&0b00001111);
            Alignment::nuc_set_as(&(nodes[z]->get_seq()->at(pos)), C[zind][i]); //preserve missing status
        }
    }
    if (profiling) Profiler::prof.stop(__func__);
    return (res);
}

double Tree::likelihood() {
    if (profiling) Profiler::prof.start(__func__);
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

    for (size_t pos = 0; pos < ali->get_npat(); pos++) {
        //Pass from leaves to root (included)
        for (size_t z = 0; z < nodes.size(); z++) {
            size_t zind=nodes[z]->get_index();
            for (int a = 0; a < 4; a++) {
                if (nodes[z]->get_left() == nullptr) {
                    if (Alignment::nuc_is_unsure((*nodes[z]->get_seq())[pos])) {
                        L[zind][a] = 1;
                    } else {
                        if (a == __builtin_ctz((*nodes[z]->get_seq())[pos]))
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
        res+=log(0.25*(L[root][0]+L[root][1]+L[root][2]+L[root][3]))*(*ali->get_weights())[pos];
    }
    if (profiling) Profiler::prof.stop(__func__);
    return res;
}

Tree* Tree::noise_tree(Tree * t,double eps) {
    size_t n=t->get_n();
    auto d=t->distmat();
    auto names=new vector<string>();
    for (size_t i=0;i<n;i++) names->push_back(t->get_node(i)->get_name());
    std::normal_distribution rnorm(0.0,eps);
    for (size_t i = 0; i < n; i++) for (size_t j = 0; j < i; j++) {
        (*d)[i][j] = abs((*d)[i][j]+rnorm(rng));
        (*d)[j][i] = (*d)[i][j];
    }
    auto t2=nj(d,names);
    delete names;
    delete d;
    return t2;
}

Tree * Tree::copy_tree(Tree * t) {
    if (profiling) Profiler::prof.start(__func__);
    Tree * t2=new Tree();
    size_t n=t->n;
    t2->n=n;
    t2->ali=t->ali;
    t2->nodes=vector<Node *>(n+n-1);
    for (size_t i = 0; i < n+n-1; i++) {
        auto node = new Node();
        node->set_name(t->get_node(i)->get_name());
        node->set_dist_father(t->get_node(i)->get_dist_father());
        node->set_index(i);
        node->set_seq(new vector<uint8_t>(*t->get_node(i)->get_seq()));
        t2->nodes[i]=node;
    }
    Node * assign;
    for (size_t i = 0; i < n+n-1; i++) {
        if (i==n+n-2) assign=nullptr; else assign=t2->get_node(t->get_node(i)->get_father()->get_index());
        t2->get_node(i)->set_father(assign);
        if (i<n) assign=nullptr; else assign=t2->get_node(t->get_node(i)->get_left()->get_index());
        t2->get_node(i)->set_left(assign);
        if (i<n) assign=nullptr; else assign=t2->get_node(t->get_node(i)->get_right()->get_index());
        t2->get_node(i)->set_right(assign);
    }
    if (profiling) Profiler::prof.stop(__func__);
    return t2;
}


void Tree::rSPR() {
    if (profiling) Profiler::prof.start(__func__);
    //First select which node is moving, which can be any node other than the root or null node
    std::uniform_int_distribution unif(0,2*static_cast<int>(n)-4);
    size_t moving=unif(rng);
    //Next select the new brother which can be any node other than root, null node, moving or its descendents
    size_t newbro;
    bool checked=false;
    while (checked==false) {
        newbro=unif(rng);
        Node*cur=nodes[newbro];
        checked=true;
        while (cur!=nullptr) {
            if (cur->get_index()==moving) {checked=false;break;}
            cur=cur->get_father();
        }
        if (newbro==nodes[moving]->get_father()->get_index()) checked=false;
    }
    //Pruning
    size_t oldf=nodes[moving]->get_father()->get_index();
    size_t oldbro=nodes[oldf]->get_left()->get_index();
    if (oldbro==moving) oldbro=nodes[oldf]->get_right()->get_index();
    double d=nodes[oldbro]->get_dist_father();
    if (nodes[oldf]->get_father()!=nullptr) {
        d+=nodes[oldf]->get_dist_father();
        int oldgf=nodes[oldf]->get_father()->get_index();
        if (nodes[oldgf]->get_left()->get_index()!=oldf) nodes[oldgf]->set_right(nodes[oldbro]); else nodes[oldgf]->set_left(nodes[oldbro]);
    }
    nodes[oldbro]->set_father(nodes[oldf]->get_father());
    nodes[oldbro]->set_dist_father(d);
    //nodes[moving]->set_dist_father(optim);//currently not optimising new branch length
    //Regrafting
    if (nodes[newbro]->get_father()!=nullptr) {
        int newbrof=nodes[newbro]->get_father()->get_index();
        if (nodes[newbrof]->get_left()->get_index()!=newbro) nodes[newbrof]->set_right(nodes[oldf]); else nodes[newbrof]->set_left(nodes[oldf]);
    }
    nodes[oldf]->set_father(nodes[newbro]->get_father());
    double m=nodes[newbro]->get_dist_father();
    std::uniform_real_distribution unifm(0.0,m);
    double x=unifm(rng);
    nodes[oldf]->set_dist_father(x);
    nodes[oldf]->set_left(nodes[moving]);
    nodes[oldf]->set_right(nodes[newbro]);
    nodes[newbro]->set_father(nodes[oldf]);
    nodes[newbro]->set_dist_father(m-x);
    //Update node matrix
    reorder();
    if (profiling) Profiler::prof.stop(__func__);
}
