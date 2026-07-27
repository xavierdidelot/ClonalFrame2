//
// Created by Didelot, Xavier on 29/09/2025.

#include "Tree.h"
#include "Alignment.h"
#include "Param.h"

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

string Tree::to_string_internal() {
    for (size_t i=n;i<nodes.size();i++) nodes[i]->set_name("NODE"+std::to_string(i));
    Node * root=nodes[2*n-2];
    return "("+root->get_left()->to_string() +","+root->get_right()->to_string() + ")"+root->get_name()+";";
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
        if (nodes[i]->get_dist_father()<0.0) {
            std::println(std::cerr,"Node with negative distance");
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
    std::uniform_int_distribution unif_l(0,static_cast<int>(l)-(l>0?1:0));
    for (size_t i = 0; i < l; i++) {
        nodes[2*n-2]->get_seq()->at(i) = Alignment::random_nuc();
    }
    int m=0;
    for (int z=2*(int)n-3;z>=0;z--) {
        vector<uint8_t>*above=nodes[z]->get_father()->get_seq();
        nodes[z]->set_seq(new vector<uint8_t>(*above));
        int m2=0;
        if (nodes[z]->get_dist_father()*l>0.0) {
                std::poisson_distribution pois(nodes[z]->get_dist_father()*l);
                m2=pois(rng);
        }
        for (int i=0;i<m2;i++) {
            int p=unif_l(rng);
            nodes[z]->get_seq()->at(p)=Alignment::random_nuc_diff(nodes[z]->get_seq()->at(p));
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
    std::println(out,"\t\t\"recombination_events\": [");
    bool first_write=true;
    for (int z=2*(int)n-3;z>=0;z--) {
        vector<uint8_t>*above=nodes[z]->get_father()->get_seq();
        nodes[z]->set_seq(new vector<uint8_t>(*above));
        double len=nodes[z]->get_dist_father();
        double a=1.0-exp(-rho*len);
        double b;
        if (a==0.0) b=0.5; else b=1.0-exp(-rho*len/expm1(rho*len*delta));//Using M/M/inf queue
        bool rec=runif(rng)<a/(a+b);
        size_t curpos=0;
        while (curpos<l) {
            size_t upto=l;
            int m2=0;
            if (len*l>0.0) {
                if ((rec?b:a)>0.0) {
                    std::geometric_distribution geo(rec?b:a);
                    upto=curpos+geo(rng)+1;
                    if (upto>l) upto=l;
                }
                std::poisson_distribution pois((rec?nu:len)*(upto-curpos));
                m2=pois(rng);
            }
            for (int i=0;i<m2;i++) {
                int p=curpos+(rng()%(upto-curpos));
                nodes[z]->get_seq()->at(p)=Alignment::random_nuc_diff(nodes[z]->get_seq()->at(p));
            }
            m+=m2;
            if (rec) {
                if (!first_write) {
                    std::println(out, ",");
                } else first_write=false;
                std::print(out, "\t\t{{\"branch\": \"{}\", \"beg\": {}, \"end\": {} }}",nodes[z]->get_name(),curpos,upto-1);
            }
            rec=!rec;
            curpos=upto;
        }
    }
    std::println(out,"]");
    simulate_pattern_analysis();
    return(m);
}

void Tree::simulate_pattern_analysis() {
    size_t l=nodes[0]->get_seq()->size();
    Alignment* a  = Alignment::make_from_tree(this,true);
    Alignment* a2 = Alignment::make_from_tree(this,false);
    a2->pattern_analysis();
    vector<int> anti(a2->get_npat(),-1);
    for (int j=static_cast<int>(l)-1;j>=0;j--) anti[a2->get_patindex()->at(j)]=j;
    for (size_t i=0;i<nodes.size();i++) {
        delete nodes[i]->get_seq();
        nodes[i]->set_seq(new vector<uint8_t>(a2->get_npat()));
        for (size_t j=0;j<a2->get_npat();j++) if (anti[j]>=0) nodes[i]->get_seq()->at(j)=a->get_data(i,anti[j]);
    }
    delete(a);
    ali=a2;
}

Tree* Tree::random_tree(size_t n,double b_mean, double b_sd) {
    auto t = new Tree();
    t->n=n;
    t->nodes=vector<Node *>(n+n-1);
    for (size_t i = 0; i < 2*n-1; i++) {
        t->nodes[i]=new Node();
        if (i<n) t->nodes[i]->set_name("seq" + std::to_string(i + 1));
    }
    auto to_merge = vector<Node *>(n);
    std::gamma_distribution rgamma(b_mean*b_mean/b_sd/b_sd,b_sd*b_sd/b_mean);
    for (size_t i = 0; i < n; i++) {
        to_merge[i]=t->nodes[i];
        to_merge[i]->set_dist_father(rgamma(rng));
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
        if (to_merge.size()>3) t->nodes[m]->set_dist_father(rgamma(rng));
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
    for (size_t i=0;i<t->nodes.size()-2;i++)
        if (std::isnan(t->nodes[i]->get_dist_father()) || t->nodes[i]->get_dist_father()<=1e-6) t->nodes[i]->set_dist_father(1e-6);
    return t;
}

void Tree::set_sequences(Alignment *ali) {
    for (size_t i=0;i<nodes.size();i++) {
        delete nodes[i]->get_seq();
        nodes[i]->set_seq(nullptr);
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
        for (size_t i=0;i<nodes.size()-2;i++) {
            res+=nodes[i]->get_dist_father();
        }
    return res;
}

double Tree::likelihood_anc() {
    if (profiling) Profiler::prof.start(__func__);
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
    if (profiling) Profiler::prof.stop(__func__);
    return ll;
}

double Tree::likelihood_anc_rec(double rho, double delta, double nu) {
    size_t L=ali->get_l();
    if (L==0) return 0.0;
    if (rho==0 || delta==0) return likelihood_anc();
    if (profiling) Profiler::prof.start(__func__);
    double ll=0.0;
    double ex_nu=exp(-4.0*nu/3.0);
    double pid_nu=0.25+0.75*ex_nu;
    double pnid_nu=0.25-0.25*ex_nu;
    vector<size_t> * refsites=ali->get_refsites();
    vector<size_t> * refsites_diff=ali->get_refsites_diff();
    vector<size_t> * patindex=ali->get_patindex();
    size_t L2=refsites->size();
    #pragma omp parallel
    {
        vector gammas(maxskip, vector(4, 0.0));
        #pragma omp for reduction(+:ll) schedule(dynamic)
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
            if (a==0.0) b=0.5; else b=1.0-exp(-rho*len/expm1(rho*len*delta));//Using M/M/inf queue
            double alpha1=b/(a+b);
            double alpha2=a/(a+b);
            vector gamma={1.0-a, a, b, 1.0-b};
            gammas[0]=gamma;
            for (size_t k=1;k<maxskip;k++) {
                gammas[k][0]=gammas[k-1][0]*pid*gamma[0]+gammas[k-1][1]*pid_nu*gamma[2];
                gammas[k][1]=gammas[k-1][0]*pid*gamma[1]+gammas[k-1][1]*pid_nu*gamma[3];
                gammas[k][2]=gammas[k-1][2]*pid*gamma[0]+gammas[k-1][3]*pid_nu*gamma[2];
                gammas[k][3]=gammas[k-1][2]*pid*gamma[1]+gammas[k-1][3]*pid_nu*gamma[3];
            }
            double lscale=0.0;
            size_t j=0;
            while (j<L2) {
                size_t skipping=(*refsites_diff)[j];
                size_t pat=(*patindex)[(*refsites)[j]];
                bool is_same=Alignment::nuc_is_currently_same((*seq1)[pat],(*seq2)[pat]);
                while (is_same) {
                    if (j==L2-1 || skipping+(*refsites_diff)[j+1]>maxskip-2) break;
                    j++;
                    skipping+=(*refsites_diff)[j]+1;
                    pat=(*patindex)[(*refsites)[j]];
                    is_same=Alignment::nuc_is_currently_same((*seq1)[pat],(*seq2)[pat]);
                }
                double p1,p2;
                if (is_same) {p1=pid;p2=pid_nu;} else {p1=pnid;p2=pnid_nu;}
                double newalpha1=(alpha1*gammas[skipping][0]+alpha2*gammas[skipping][2])*p1;
                alpha2=          (alpha1*gammas[skipping][1]+alpha2*gammas[skipping][3])*p2;
                alpha1=newalpha1;
                double scale=alpha1+alpha2;
                alpha1=alpha1/scale;
                alpha2=alpha2/scale;
                lscale+=log(scale);
                j++;
            }
            ll+=lscale;
        }
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
    todo[0]->set_dist_father(0.0);
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
    return d;
}

double Tree::jointml() {
    if (ali->get_l()==0) return 0.0;
    if (profiling) Profiler::prof.start(__func__);
    double res = 0.0;
    const double log025 = log(0.25);

    vector<double> pid(nodes.size(),0);
    vector<double> pnid(nodes.size(),0);
    for (size_t z = 0; z < n*2-1; z++) {
        size_t zind=nodes[z]->get_index();
        double ex = exp(-4.0*nodes[z]->get_dist_father()/3.0);
        pid[zind] = log(0.25 + 0.75 * ex);
        pnid[zind] = log(0.25 - 0.25 * ex);
    }

    #pragma omp parallel
    {
        vector<vector<double> > L(nodes.size(), vector<double>(4, 0.0));
        vector<vector<uint8_t> > C(nodes.size(), vector<uint8_t>(4, 0));
        double tmp[4];
        int z, w;
        #pragma omp for reduction(+:res) schedule(static)
        for (size_t pos = 0; pos < ali->get_npat(); pos++) {
            //First pass for leaves
            for (z = 0; z < static_cast<int>(n); z++) {
                size_t zind = nodes[z]->get_index();
                uint8_t nuc = nodes[z]->get_seq()->at(pos);
                int nucint = __builtin_ctz(nuc);
                bool nuc_is_unsure = Alignment::nuc_is_unsure(nuc);
                for (int i = 0; i < 4; i++) {
                    if (nuc_is_unsure) {
                        L[zind][i] = pid[zind];
                        C[zind][i] = 1 << i;
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
            for (z = static_cast<int>(n); z < static_cast<int>(n) * 2 - 1; z++) {
                size_t zind = nodes[z]->get_index();
                size_t left = nodes[z]->get_left()->get_index();
                size_t right = nodes[z]->get_right()->get_index();
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
            z = 2 * n - 2;
            size_t left = nodes[z]->get_left()->get_index();
            size_t right = nodes[z]->get_right()->get_index();
            for (size_t j = 0; j < 4; j++) {
                L[z][j] = L[left][j] + L[right][j] + log025;
            }
            w = 0;
            for (size_t j = 1; j < 4; j++) { if (L[z][j] > L[z][w]) w = j; }
            res += L[z][w] * ali->get_weights()->at(pos);
            Alignment::nuc_set_as(&(nodes[z]->get_seq()->at(pos)), 1 << w); //preserve missing status

            //Second pass from root (excluded) to leaves
            for (z = 2 * n - 3; z >= 0; z--) {
                size_t zind = nodes[z]->get_index();
                size_t i = __builtin_ctz(nodes[z]->get_father()->get_seq()->at(pos) & 0b00001111);
                Alignment::nuc_set_as(&(nodes[z]->get_seq()->at(pos)), C[zind][i]); //preserve missing status
            }
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

double Tree::prior(Param * p) {
    auto prior_params=*(p->get_prior_params());
    if (prior_params[6]==prior_params[7]) {
        //Special Exponential case
        return -sum_bra_len()/prior_params[6]-log(prior_params[6])*(2.0*n-3.0);
    }
    double shape=prior_params[6]*prior_params[6]/prior_params[7]/prior_params[7];
    double scale=prior_params[7]*prior_params[7]/prior_params[6];
    double res=0;
    for (size_t i=0;i<nodes.size()-2;i++) {
        res+=dgamma_shape_scale(nodes[i]->get_dist_father(),shape,scale);
    }
    return res;
}

Tree* Tree::noise_tree(Tree * t,double eps) {
    size_t n=t->get_n();
    auto d=t->distmat();
    auto names=new vector<string>();
    for (size_t i=0;i<n;i++) names->push_back(t->get_node(i)->get_name());
    std::normal_distribution rnorm(0.0,eps);
    for (size_t i = 0; i < n; i++) for (size_t j = 0; j < i; j++) {
        (*d)[i][j] = std::abs((*d)[i][j]+rnorm(rng));
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

double Tree::select_rSPR(size_t moving, Param * p, size_t * newbro, double * m, double * f, double * b, double * d) {
    auto prior_params=*(p->get_prior_params());
    size_t max=2*static_cast<int>(n)-3;
    vector probs=vector<double>(max+1,std::log(0.0));
    //Select the new brother which can be any node other than root, null node, moving, its father or its descendents
    for (size_t newb=0;newb<max;newb++) {
        Node*cur=nodes[newb];
        bool checked=true;
        while (cur!=nullptr) {
            if (cur->get_index()==moving) {checked=false;break;}
            cur=cur->get_father();
        }
        if (newb==nodes[moving]->get_father()->get_index()) checked=false;
        if (checked) {
            int diffs=count_diffs(nodes[moving],nodes[newb]);
            probs[newb]=ali->get_l()==0?0:std::log(ali->get_l()*prior_params[6]/(1.0+ali->get_l()*prior_params[6]))*diffs;
        }
    }
    double probs_max = *std::max_element(probs.begin(), probs.end());
    for (double &x : probs) x = std::exp(x-probs_max);
    size_t oldbro=nodes[moving]->get_brother()->get_index();
    double sum_minus_old = std::accumulate(probs.begin(), probs.end(), 0.0)-probs[oldbro];
    if (sum_minus_old>0.0) {
        for (double &x : probs) x=x/sum_minus_old/2.0;
        probs[oldbro]=0.5;
    } else probs[oldbro]=1.0;
    std::discrete_distribution<int> dist(probs.begin(), probs.end());
    *newbro=dist(rng);
    //*newbro=oldbro;//Use only for testing: forces the regrafting to be on the same branch as before
    double q=log(probs[*newbro]);//This is the probability of having selected newbro

    double shape=1,scale=1;
    double d1,d2,d3,dm,df,db;
    Node * newfather=nodes[*newbro]->get_father();
    if (*newbro==oldbro) newfather=newfather->get_father();
    d1=count_diffs(nodes[moving],newfather);
    d2=count_diffs(nodes[moving],nodes[*newbro]);
    d3=count_diffs(nodes[*newbro],newfather);
    dm=(d1+d2-d3)/2.0;
    df=(d1+d3-d2)/2.0;
    db=(d2+d3-d1)/2.0;
    //if (verbose) vprint("d1={:.2e} d2={:.2e} d3={:.2e} dm={:.2e} df={:.2e} db={:.2e}",d1,d2,d3,dm,df,db);
    std::gamma_distribution<double> rgamma(shape, scale);
    determine_shape_scale(dm,p,&shape,&scale);
    rgamma.param(std::gamma_distribution<double>::param_type(shape,scale));
    *m=rgamma(rng);
    q+=dgamma_shape_scale(*m,shape,scale);//This is the dgamma probability of the new distance m

    determine_shape_scale(df,p,&shape,&scale);
    rgamma.param(std::gamma_distribution<double>::param_type(shape,scale));
    *f=rgamma(rng);
    q+=dgamma_shape_scale(*f,shape,scale);//This is the dgamma probability of the new distance f

    determine_shape_scale(db,p,&shape,&scale);
    rgamma.param(std::gamma_distribution<double>::param_type(shape,scale));
    *b=rgamma(rng);
    q+=dgamma_shape_scale(*b,shape,scale);//This is the dgamma probability of the new distance b

    if (*newbro!=nodes[moving]->get_brother()->get_index()) {
        int diffs=count_diffs(nodes[moving]->get_brother(),nodes[moving]->get_father()->get_father());
        determine_shape_scale(diffs,p,&shape,&scale);
        rgamma.param(std::gamma_distribution<double>::param_type(shape,scale));
        *d=rgamma(rng);
        q+=dgamma_shape_scale(*d,shape,scale);//This is the dgamma probability of the new distance d
    } else *d=-1.0;
    return(q);
}

double Tree::rSPR(Param * p,string * name_moving) {
    if (n==3) return 0;
    if (profiling) Profiler::prof.start(__func__);
    //First select which node is moving, which can be any node other than the root, the null node and their three children
    std::uniform_int_distribution unif(0,2*static_cast<int>(n)-4);
    size_t moving,oldf,oldbro;
    bool checked=false;
    while (checked==false) {
        moving=unif(rng);
        oldf=nodes[moving]->get_father()->get_index();
        oldbro=nodes[moving]->get_brother()->get_index();
        if (oldf<=2*n-4 && oldbro<=2*n-4) checked=true;
    }
    //Next select where and how to regraft
    double m,f,b,d;
    size_t newbro;
    double q=select_rSPR(moving, p,&newbro, &m, &f, &b, &d);
    //Pruning
    Node * oldgf=nodes[oldf]->get_father();
    if (oldgf->get_left()->get_index()==oldf) oldgf->set_left(nodes[oldbro]); else oldgf->set_right(nodes[oldbro]);
    nodes[oldbro]->set_father(nodes[oldf]->get_father());
    //if (verbose) vprint("bro    {:.2e} {:.2e} {:.2e}",nodes[oldbro]->get_dist_father(),b,nodes[oldbro]->get_dist_father()-b);//because d=-1 so anticipating next change
    nodes[oldbro]->set_dist_father(d);
    //if (verbose) vprint("moving {:.2e} {:.2e} {:.2e}",nodes[moving]->get_dist_father(),m,nodes[moving]->get_dist_father()-m);
    nodes[moving]->set_dist_father(m);
    //Regrafting
    Node * newbrof=nodes[newbro]->get_father();
    if (newbrof->get_left()->get_index()!=newbro) newbrof->set_right(nodes[oldf]); else newbrof->set_left(nodes[oldf]);
    nodes[oldf]->set_father(nodes[newbro]->get_father());
    //if (verbose) vprint("father {:.2e} {:.2e} {:.2e}",nodes[oldf]->get_dist_father(),f,nodes[oldf]->get_dist_father()-f);
    nodes[oldf]->set_dist_father(f);
    nodes[oldf]->set_left(nodes[moving]);
    nodes[oldf]->set_right(nodes[newbro]);
    nodes[newbro]->set_father(nodes[oldf]);
    nodes[newbro]->set_dist_father(b);
    //Update node matrix
    *name_moving=nodes[moving]->get_name();
    reorder();
    if (profiling) Profiler::prof.stop(__func__);
    return(q);
}

double Tree::lq_rSPR(Tree *tree2, Param * p,string name_moving) {
    auto prior_params=*(p->get_prior_params());
    if (n==3) return 0;
    if (profiling) Profiler::prof.start(__func__);
    size_t moving=0;
    while (this->nodes[moving]->get_name()!=name_moving) moving++;
    size_t moving2=0;
    while (tree2->nodes[moving2]->get_name()!=name_moving) moving2++;
    string name_newbro=tree2->nodes[moving2]->get_brother()->get_name();
    size_t newbro=0;
    while (this->nodes[newbro]->get_name()!=name_newbro) newbro++;
    string name_oldbro=this->nodes[moving]->get_brother()->get_name();
    size_t oldbro2=0;
    while (tree2->nodes[oldbro2]->get_name()!=name_oldbro) oldbro2++;
    size_t oldbro=0;
    while (this->nodes[oldbro]->get_name()!=name_oldbro) oldbro++;
    double m=tree2->nodes[moving2]->get_dist_father();
    double f=tree2->nodes[moving2]->get_father()->get_dist_father();
    double b=tree2->nodes[moving2]->get_brother()->get_dist_father();
    double d=tree2->nodes[oldbro2]->get_dist_father();

    double q=0.0;
    size_t max=2*static_cast<int>(n)-3;
    vector probs=vector<double>(max+1,std::log(0.0));
    for (size_t newb=0;newb<max;newb++) {
        Node*cur=nodes[newb];
        bool checked=true;
        while (cur!=nullptr) {
            if (cur->get_index()==moving) {checked=false;break;}
            cur=cur->get_father();
        }
        if (newb==nodes[moving]->get_father()->get_index()) checked=false;
        if (checked) {
            int diffs=count_diffs(nodes[moving],nodes[newb]);
            probs[newb]=ali->get_l()==0?0:std::log(ali->get_l()*prior_params[6]/(1.0+ali->get_l()*prior_params[6]))*diffs;
        }
    }
    double probs_max = *std::max_element(probs.begin(), probs.end());
    for (double &x : probs) x = std::exp(x-probs_max);
    double sum_minus_old = std::accumulate(probs.begin(), probs.end(), 0.0)-probs[oldbro];
    if (sum_minus_old>0.0) {
        for (double &x : probs) x=x/sum_minus_old/2.0;
        probs[oldbro]=0.5;
    } else probs[oldbro]=1.0;
    q+=log(probs[newbro]);//This is the probability of having selected newbro
    double shape,scale;
    double d1,d2,d3,dm,df,db;
    Node * newfather=nodes[newbro]->get_father();
    if (newbro==oldbro) newfather=newfather->get_father();
    d1=count_diffs(nodes[moving],newfather);
    d2=count_diffs(nodes[moving],nodes[newbro]);
    d3=count_diffs(nodes[newbro],newfather);
    dm=(d1+d2-d3)/2.0;
    df=(d1+d3-d2)/2.0;
    db=(d2+d3-d1)/2.0;
    determine_shape_scale(dm,p,&shape,&scale);
    q+=dgamma_shape_scale(m,shape,scale);//This is the dgamma probability of the new distance m
    determine_shape_scale(df,p,&shape,&scale);
    q+=dgamma_shape_scale(f,shape,scale);//This is the dgamma probability of the new distance f
    determine_shape_scale(db,p,&shape,&scale);
    q+=dgamma_shape_scale(b,shape,scale);//This is the dgamma probability of the new distance b
    if (nodes[newbro]->get_name()!=nodes[oldbro]->get_name()) {
        int diffs=count_diffs(nodes[oldbro],nodes[oldbro]->get_father()->get_father());
        determine_shape_scale(diffs,p,&shape,&scale);
        q+=dgamma_shape_scale(d,shape,scale);//This is the dgamma probability of the new distance d
    }
    if (profiling) Profiler::prof.stop(__func__);
    return q;
}

void Tree::determine_shape_scale(double diffs, Param * p, double *shape, double *scale) {
    auto prior_params=*(p->get_prior_params());
    double b_shape=prior_params[6]*prior_params[6]/prior_params[7]/prior_params[7];
    double b_scale=prior_params[7]*prior_params[7]/prior_params[6];
    double rho=p->get_rho();
    double delta=p->get_delta();
    double nu=p->get_nu();
    nu=0.75-0.75*std::exp(-4.0*nu/3.0);
    double L=ali->get_l();
    diffs=diffs/(1.0+rho*delta*nu);//Corrected diffs for effect of recombination
    double d=0.0;
    if (L>0) d=-0.75*std::log(1.0-4.0/3.0*diffs/L);
    diffs=d*L;
    double lambda=2;
    *shape=b_shape+diffs/lambda;
    *scale=b_scale*lambda/(L*b_scale+lambda);
}

int Tree::count_diffs(Node *n1, Node *n2) {
    int diffs=0;
    for (size_t pos=0;pos<ali->get_npat();pos++) {
        if (!Alignment::nuc_is_currently_same(n1->get_seq()->at(pos),n2->get_seq()->at(pos))) diffs+=ali->get_weights()->at(pos);
    }
    return(diffs);
}

void Tree::ancestral_to_json_file() {
    auto patindex=ali->get_patindex();
    auto npat=ali->get_npat();
    auto l=ali->get_l();
    std::println(out, "\t\"ancestral_sequences\": {{");
    std::print(out, "\t\t\"patindex\": [");
    if (l>0) std::print(out, "{}",patindex->at(0));
    for (size_t i=1;i<l;i++) {
        std::print(out, ",{}",patindex->at(i));
    }
    std::println(out, "],");
    std::println(out, "\t\t\"patmatrix\": [");
    for (size_t i=0;i<2*n-2;i++) {
        std::print(out, "\t\t{{\"branch\": \"{}\", \"data\": \"",nodes[i]->get_name());
        for (size_t j=0;j<npat;j++) std::print(out,"{}",Alignment::decode(0b00001111&nodes[i]->get_seq()->at(j)));
        if (i<2*n-3) std::println(out, "\"}},"); else std::println(out, "\"}}");
    }
    std::println(out, "\t\t]");
    std::println(out, "\t}},");
}
