#include "../src/Alignment.h"
#include <cassert>

using std::string;

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
    verbose=true;
    std::ofstream out("tmp.fasta");
    out << fasta;
    out.close();
    auto a = new Alignment("tmp.fasta");
    remove("tmp.fasta");
    a->write("tmp2.fasta");
    remove("tmp2.fasta");
    size_t n=a->get_n();assert(n==3);
    size_t l=a->get_l();assert(l==8);
    size_t npat=a->get_npat();assert(npat==5);
    uint8_t c=a->get_data(1,1);assert(c==2);
    string s=a->get_name(2);assert(s=="seq3");
    delete(a);
}