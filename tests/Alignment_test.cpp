#include "../src/Alignment.h"
#include <cassert>

static string fasta=
    R"(
>seq1
ACGTACGTACGT
> seq2
acgtACGTACGA
>  seq3
ACGTACGANNNN
)";


int main() {
    verbose=true;
    std::ofstream out("tmp.fasta");
    myprintln(out,"{}",fasta);
    out.close();
    auto a = Alignment::read_fasta("tmp.fasta");
    remove("tmp.fasta");
    a->write("tmp2.fasta",true);
    remove("tmp2.fasta");
    size_t n=a->get_n();assert(n==3);
    size_t l=a->get_l();assert(l==12);
    size_t npat=a->get_npat();assert(npat==5+1);
    uint8_t c=a->get_data(1,1);assert(c==2);
    string s=a->get_name(2);assert(s=="seq3");
    delete(a);
}