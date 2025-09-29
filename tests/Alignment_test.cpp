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
    assert(a->get_n()==3);
    assert(a->get_l()==8);
    assert(a->get_data(1,1)=='C');
    assert(a->get_name(2)=="seq3");
    delete(a);
}