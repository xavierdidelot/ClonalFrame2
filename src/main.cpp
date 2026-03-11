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
    -s NUM      Use given seed to initiate random number generator
    -S          Generate a simulated dataset and run on it
    -v          Verbose mode
)";

int main(int argc, const char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    auto *orphan = new Node();orphan->set_dist_father(1.0); //just to test leak detection is working

    int burn = 10000, nbit = 10000, thin = 10;
    int i=1;
    bool simu=false;
    vector<string> remainder;
    while (i < argc) {
        string arg = argv[i];
        if (arg[0]=='-') {
            switch (arg[1]) {
                case('x'):burn =  std::stoi(argv[++i]) ;break;
                case('y'):nbit =  std::stoi(argv[++i]) ;break;
                case('z'):thin =  std::stoi(argv[++i]) ;break;
                case('s'):setSeed(std::stoi(argv[++i]));break;
                case('S'):simu = true;break;
                case('v'):verbose = true;break;
                default:cout << "Ignoring incorrect option -" << arg[1] << endl;
            }
        } else remainder.push_back(arg);
        i++;
    }

    if (remainder.size() != 2 && simu==false) {
        cout << "This is ClonalFrame2, version " << VERSION_STRING << endl;
        cout << "Need two non-option arguments but found " << remainder.size() << endl;
        cout << help << endl;
    } else {
        if (verbose) {
            cout << "This is ClonalFrame2, version " << VERSION_STRING << endl;
            cout << "Command used: ";
            for (i = 0; i < argc; i++) cout << argv[i] << " ";
            cout << endl;
        }

        Alignment * a;
        string outfile;
        if (simu) {
            if (verbose) cout << "Generating true tree..." << endl;
            auto t0 = Tree::random_tree(100,1e-4);
            if (verbose) cout << "Generating simulated alignment..." << endl;
            t0->simulate_sequences(1e6);
            a = Alignment::make_from_tree(t0,false,false);
            delete t0;
            outfile="";
        } else {
            const auto inputfile = string(remainder[0]);
            outfile = string(remainder[1]);
            a = Alignment::read_fasta(inputfile);
        }
        if (verbose) cout << "Building initial tree..." << endl;
        auto t = Tree::nj(a);
        t->set_sequences(a);
        if (verbose) cout << "Performing initial ancestral reconstruction..." << endl;
        t->jointml();
        auto p = new Param(t);
        p->mcmc(burn, nbit, thin, outfile);
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