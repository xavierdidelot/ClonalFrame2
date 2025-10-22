#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

using namespace std;

int main() {
    auto t = new Tree(10);
    size_t n=t->get_n();assert(n==10);
    bool b=t->check_structure();assert(b);
    t->randomise(1.0);
    b=t->check_structure();assert(b);
    t->write("tmp.nwk");
    auto t2=new Tree("tmp.nwk");
    b=t2->check_structure();assert(b);
    t2->write("t2.nwk");
    delete t2;
    t->add_mutations(10000,10);
    auto a=new Alignment(t);
    a->miss_at_random(0.1);
    a->write("tmp.fasta");
    int truth=t->count_subs();
    cout<<"truth="<<truth<<endl;
    t->nj(a);
    b=t->check_structure();assert(b);
    t->write("tmp2.nwk");
    t->set_seqs(a);
    int s1=t->ancestral_reconstruction();
    int s2=t->count_subs();
    cout<<s1<<" "<<s2<<endl;
    assert(s1==s2);
    delete(t);
    delete(a);
}