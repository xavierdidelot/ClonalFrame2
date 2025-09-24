#include "../src/Alignment.h"
#include <cassert>

static string fasta=
    R"(
>1
ACGT
>2
ACGT
>3
ACGT
)";


int main() {
    std::ofstream out("tmp.fasta");
    out << fasta;
    out.close();
    auto *a = new Alignment("tmp.fasta");
    assert(a->getN()==3);
    assert(a->getL()==4);
    delete(a);
    remove("tmp.fasta");
}