#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

//This testing main function reads the alignment in tmp.fasta and count the number of unique patterns
int main() {
    auto a = Alignment::read_fasta("tmp.fasta");
    //Count number of patterns actually used rather than rely on a->get_npat() which always includes non-polymorphic patterns
    vector<size_t> * pat=a->get_patindex();
    auto tmp = *pat;
    std::sort(tmp.begin(), tmp.end());
    size_t count = std::unique(tmp.begin(), tmp.end()) - tmp.begin();
    myprintln("{}",count);
    delete(a);
}