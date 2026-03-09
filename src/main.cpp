#include <iostream>
#include <chrono>

#include "Alignment.h"
#include "Tree.h"
#include "Param.h"
#include "RNGGlobal.h"
#include "version.h"

using std::string,std::vector,std::cerr,std::cout,std::endl;

static string help =
        R"(Usage: ClonalFrame2 [OPTIONS] inputfile outputfile

    Options:
    -x NUM      Sets the number of burn-in iterations (default is 10000)
    -y NUM      Sets the number of iterations after burn-in (default is 10000)
    -z NUM      Sets the number of iterations between samples (default is 10)
    -s NUM      Use given seed to initial random number generator
    -v          Verbose mode
)";

int main(int argc, const char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    auto *orphan = new Node();orphan->set_dist_father(1.0); //just to test leak detection is working

    int burn = 10000, nbit = 10000, thin = 10;
    int i=1;
    vector<string> remainder;
    while (i < argc) {
        string arg = argv[i];
        if (arg[0]=='-') {
            switch (arg[1]) {
                case('x'):burn =  std::stoi(argv[++i]) ;break;
                case('y'):nbit =  std::stoi(argv[++i]) ;break;
                case('z'):thin =  std::stoi(argv[++i]) ;break;
                case('s'):setSeed(std::stoi(argv[++i]));break;
                case('v'):verbose = true;break;
                default:cout << "Ignoring incorrect option -" << arg[1] << endl;
            }
        } else remainder.push_back(arg);
        i++;
    }

    if (remainder.size() != 2) {
        cout << "Need two non-option arguments but found " << remainder.size() << endl;
        cout << help << endl;
    } else {
        if (verbose) {
            cout << "This is ClonalFrame2, version " << VERSION_STRING << endl;
            cout << "Command used: ";
            for (i = 0; i < argc; i++) cout << argv[i] << " ";
            cout << endl;
        }
        const auto inputfile = string(remainder[0]);
        const auto outfile = string(remainder[1]);
        auto a = Alignment::read_fasta(inputfile);
        if (verbose) cout << "Building initial tree..." << endl;
        auto t = Tree::nj(a);
        t->set_sequences(a);
        if (verbose) cout << "Performing initial ancestral reconstruction..." << endl;
        t->parsimony();
        auto p = new Param(t);
        p->mcmc(burn, nbit, thin);
        t->write(outfile);
        delete(a);
        delete(t);
        delete(p);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    if (verbose) cout << "Time spent running: " << elapsed.count() << " seconds" << endl;
    if (verbose) cout << "End of ClonalFrame2 run." << endl;
    return 0;
}