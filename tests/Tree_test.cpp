#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

using namespace std;

int main() {

    //Creating a tree
    auto t = Tree::deterministic_tree(10);
    size_t n=t->get_n();assert(n==10);
    bool b=t->check_structure();assert(b);
    delete(t);

    //Random tree
    t = Tree::random_tree(10,0.01);
    n=t->get_n();assert(n==10);
    b=t->check_structure();assert(b);
    cout<<"Expected number of substitutions: "<<t->sum_bra_len()*10000<<endl;

    //Writing and re-reading a tree
    t->write("true.nwk");
    auto t2=Tree::read_newick("true.nwk");
    b=t2->check_structure();assert(b);
    delete t2;

    //Generating an alignment from a tree
    int m=t->simulate_sequences(10000);
    cout<<"Simulated number of substitutions: "<<m<<endl;
    auto a=Alignment::make_from_tree(t,false,false);
    //a->miss_at_random(0.1);
    a->write("ali.fasta",true);

    //Counting number of substitutions in the tree
    int truth=t->count_subs();
    cout<<"Counted number of substitutions: "<<truth<<endl;

    //Ancestral reconstruction
    int s1=t->parsimony();
    int s2=t->count_subs();
    int s3=t->parsimony_recursive();
    int s4=t->count_subs();
    cout<<"Parsimonious number of substitutions: "<<s1<<" "<<s2<<" "<<s3<<" "<<s4<<endl;
    assert(s1==s2);
    assert(s2==s3);
    assert(s3==s4);
    cout<<"Likelihood given parsimony ancestral sequences: "<<std::fixed << std::setprecision(6)<<t->likelihood_anc()<<endl;
    double l1=t->jointml();
    double l2=t->likelihood_anc();
    cout<<"Likelihoods given ML ancestral sequences: "<<std::fixed << std::setprecision(6)<<l1<<" "<<std::fixed << std::setprecision(6)<<l2<<endl;
    assert(round(l1*1e6)==round(l2*1e6));
    //NJ algorithm
    delete(t);
    t=Tree::nj(a);
    b=t->check_structure();assert(b);
    t->write("nj.nwk");
    t->set_sequences(a);
    delete(t);
    delete(a);
}