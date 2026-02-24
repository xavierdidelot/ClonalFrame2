#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

using namespace std;

//This testing main function reads the alignment in tmp.fasta and count the number of unique patterns
int main() {
    auto a = Alignment::read_fasta("tmp.fasta");
    cout<<a->get_npat()<<endl;
    delete(a);
}