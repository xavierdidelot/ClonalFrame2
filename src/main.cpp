#include <iostream>
#include <chrono>

#include "Alignment.h"
#include "Tree.h"
#include "RNGGlobal.h"
#include "version.h"

static string help =
        R"(Usage: ClonalFrame2 [OPTIONS] inputfile outputfile

    Options:
    -x NUM      Sets the number of burn-in iterations (default is 10000)
    -y NUM      Sets the number of iterations after burn-in (default is 10000)
    -z NUM      Sets the number of iterations between samples (default is 10)
    -s NUM      Use given seed to initial random number generator
    -v          Verbose mode
)";

bool verbose = false;

int main(int argc, char *argv[]) {
    auto start = chrono::high_resolution_clock::now();
    cout << "This is ClonalFrame2, version " << VERSION_STRING << endl;
    cout << "Command used: ";
    for (int i = 0; i < argc; i++) cout << argv[i] << " ";
    cout << endl;
    auto *orphan = new Node();orphan->set_dist_father(1.0); //just to test leak detection is working

    int burn = 10000, nbit = 10000, thin = 10;
    int i=1;
    vector<string> remainder;
    while (i < argc) {
        string arg = argv[i];
        if (arg[0]=='-') {
            char opt = arg[1];
            switch (opt) {
                case('x'):burn =  stoi(argv[++i]) ;break;
                case('y'):nbit =  stoi(argv[++i]) ;break;
                case('z'):thin =  stoi(argv[++i]) ;break;
                case('s'):setSeed(stoi(argv[++i]));break;
                case('v'):verbose = true;break;
                default:cout << "Ignoring incorrect option -" << opt << endl;
            }
        } else remainder.push_back(arg);
        i++;
    }

    if (remainder.size() != 2) {
        cout << "Need two non-option arguments but found " << remainder.size() << endl;
        cout << help << endl;
    } else {
        string inputfile = string(remainder[0]);
        string outfile = string(remainder[1]);
        swap(burn, nbit);swap(burn,thin); //just to avoid warnings of unused variables
        auto *a = new Alignment(inputfile);
        cout << "Building initial tree..." << endl;
        auto *t = new Tree(a->get_n());
        t->nj(a);
        t->set_seqs(a);
        cout << "Performing initial ancestral reconstruction..." << endl;
        t->ancestral_reconstruction();
        t->write(outfile);
        delete a;
        delete t;
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "Time spent running: " << elapsed.count() << " seconds" << endl;
    cout << "End of ClonalFrame2 run." << endl;
    return 0;
}