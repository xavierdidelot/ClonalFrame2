#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

//This testing main function reads the alignment in tmp.fasta and count the number of unique patterns
int main() {
    auto a = Alignment::read_fasta("tmp.fasta");
    std::println("{}",a->get_npat());
    delete(a);
}