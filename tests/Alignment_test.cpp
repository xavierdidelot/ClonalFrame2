#include "../src/Alignment.h"
#include <cassert>

static string fasta=
    R"(
>seq1
ACGTACGT
> seq2
acgtACGT
>  seq3
ACGTACGA
)";


int main() {
    std::ofstream out("tmp.fasta");
    out << fasta;
    out.close();
    auto *a = new Alignment("tmp.fasta");
    remove("tmp.fasta");
    a->write("tmp2.fasta");
    remove("tmp2.fasta");
    int n=a->get_n();assert(n==3);
    int l=a->get_l();assert(l==8);
    char c=a->get_data(1,1);assert(c==2);
    string s=a->get_name(2);assert(s=="seq3");
    delete(a);
}