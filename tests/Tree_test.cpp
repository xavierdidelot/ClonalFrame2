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
    t = Tree::random_tree(10,1.0);
    n=t->get_n();assert(n==10);
    b=t->check_structure();assert(b);
    cout<<"Expected number of substitutions: "<<t->sum_bra_len()*10/2<<endl;

    //Writing and re-reading a tree
    t->write("true.nwk");
    auto t2=Tree::read_newick("true.nwk");
    b=t2->check_structure();assert(b);
    delete t2;

    //Generating an alignment from a tree
    t->simulate_sequences(10000,10);
    auto a=Alignment::make_from_tree(t);
    //a->miss_at_random(0.1);
    a->write("ali.fasta");

    //Counting number of substitutions in the tree
    int truth=t->count_subs();
    cout<<"Actual number of substitutions: "<<truth<<endl;

    //Ancestral reconstruction
    int s1=t->ancestral_reconstruction();
    int s2=t->count_subs();
    cout<<"Parsimonious number of substitutions: "<<s1<<" "<<s2<<endl;
    assert(s1==s2);

    //NJ algorithm
    delete(t);
    t=Tree::nj(a);
    b=t->check_structure();assert(b);
    t->write("nj.nwk");
    t->set_sequences(a);
    delete(t);
    delete(a);
}