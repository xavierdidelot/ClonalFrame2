#include "../src/Alignment.h"
#include <cassert>

static string fasta=
    R"(
>seq1
ACGTACGT
> seq2
acgtACGT
>  seq3
ACGTACGT
)";


int main() {
    std::ofstream out("tmp.fasta");
    out << fasta;
    out.close();
    auto *a = new Alignment("tmp.fasta");
    remove("tmp.fasta");
    a->write("tmp2.fasta");
    remove("tmp2.fasta");
    assert(a->getN()==3);
    assert(a->getL()==8);
    assert(a->getData(1,1)=='C');
    assert(a->getName(2)=="seq3");
    delete(a);
}