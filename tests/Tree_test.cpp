#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

int main() {
    setSeed(0);

    //Creating a tree
    auto t = Tree::deterministic_tree(10);
    size_t n=t->get_n();assert(n==10);
    bool b=t->check_structure();assert(b);
    delete(t);

    //Random tree
    t = Tree::random_tree(10,0.01);
    n=t->get_n();assert(n==10);
    b=t->check_structure();assert(b);
    std::println("Expected number of substitutions: {}",t->sum_bra_len()*10000);

    //Writing and re-reading a tree
    t->write("true.nwk");
    auto t2=Tree::read_newick("true.nwk");
    b=t2->check_structure();assert(b);
    n=t2->get_n();assert(n==10);
    delete t2;

    //Generating an alignment from a tree
    int m=t->simulate_sequences(10000);
    std::println("Simulated number of substitutions: {}",m);
    auto a=Alignment::make_from_tree(t,false,false);
    //a->miss_at_random(0.1);
    a->write("ali.fasta",true);

    //Likelihood calculation
    std::println("Likelihood: {}",t->likelihood());

    //Counting number of substitutions in the tree
    int truth=t->count_subs();
    std::println("Counted number of substitutions: {}",truth);

    //Ancestral reconstruction
    int s1=t->parsimony();
    int s2=t->count_subs();
    int s3=t->parsimony_recursive();
    int s4=t->count_subs();
    std::println("Parsimonious number of substitutions: {} {} {} {}",s1,s2,s3,s4);
    assert(s1==s2);
    assert(s2==s3);
    assert(s3==s4);
    std::println("Likelihood given parsimony ancestral sequences: {}",t->likelihood_anc());
    double l1=t->jointml();
    double l2=t->likelihood_anc();
    std::println("Likelihoods given ML ancestral sequences: {} {}",l1,l2);
    assert(round(l1*1e6)==round(l2*1e6));

    //NJ algorithm
    delete(t);
    auto mat=a->distmat();
    t = Tree::nj(mat,a->get_names());
    delete mat;
    b=t->check_structure();assert(b);
    t->write("nj.nwk");
    t->set_sequences(a);

    //rSPR algorithm
    delete(t);
    t=Tree::random_tree(20,0.01);
    //double br=t->sum_bra_len();
    for (int i=0;i<100;i++) {
        //t->rSPR();//TODO
        b=t->check_structure();
        assert(b);
    }

    //Noisy algorithm
    auto t3=Tree::noise_tree(t,0.01);
    b=t3->check_structure();assert(b);
    t3->write("nj_noise.nwk");
    delete t3;
    delete(t);
    delete(a);

    //Distance matrix
    std::ofstream out("tmp.nwk");
    std::println(out,"{}","((seq2:1,seq1:2):1,seq0:3);");
    out.close();
    t=Tree::read_newick("tmp.nwk");
    n=t->get_n();assert(n==3);
    vector<vector<double>> * d=t->distmat();
    assert((*d)[0][0]==0);
    assert((*d)[0][1]==6);
    assert((*d)[0][2]==5);
    delete(t);
}