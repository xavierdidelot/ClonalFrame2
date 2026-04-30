#include <iostream>
#include <chrono>

#include "Alignment.h"
#include "Tree.h"
#include "Param.h"
#include "RNGGlobal.h"
#include "version.h"

static string help =
        R"(Usage: ClonalFrame2 [OPTIONS] inputfile outputfile

    Options:
    -x NUM      Sets the number of burn-in iterations (default is 10000)
    -y NUM      Total number of iterations (default is 20000)
    -z NUM      Sets the number of iterations between samples (default is 10)
    -s NUM      Use given seed to initiate random number generator
    -S          Generate a simulated dataset and run on it
    -v          Verbose mode
)";

int main(int argc, const char *argv[]) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto *orphan = new Node();orphan->set_dist_father(1.0); //just to test leak detection is working

    int burn = 10000, nbit = 20000, thin = 10;
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
                default:std::println("Ignoring incorrect option -{}",arg[1]);
            }
        } else remainder.push_back(arg);
        i++;
    }

    std::println("This is ClonalFrame2, version {}",VERSION_STRING);
    if (remainder.size() != 2 && simu==false) {
        std::println("Need two non-option arguments but found {}",remainder.size());
        std::println("{}",help);
    } else {
        if (verbose) {
            std::print("Command used: ");
            for (i = 0; i < argc; i++) std::print("{} ",argv[i]);
            std::println("");
        }

        Alignment * a;
        string outfile;
        if (simu) {
            if (verbose) std::println("Generating true tree...");
            auto t0 = Tree::random_tree(100,1e-4);
            if (verbose) std::println("Generating simulated alignment...");
            t0->simulate_sequences(1e6);
            a = Alignment::make_from_tree(t0,false);
            delete t0->get_patindex();
            delete t0->get_weights();
            delete t0;
            outfile="";
        } else {
            const auto inputfile = string(remainder[0]);
            outfile = string(remainder[1]);
            a = Alignment::read_fasta(inputfile);
        }
        if (verbose) std::println("Building initial tree...");
        auto mat=a->distmat();
        auto t = Tree::nj(mat,a->get_names());
        delete mat;
        t->set_sequences(a);
        if (verbose) std::println("Performing initial ancestral reconstruction...");
        t->jointml();
        auto p = new Param(t);
        if (verbose) std::println("Starting MCMC...");
        p->mcmc(burn, nbit, thin, outfile);
        t=p->get_tree();
        if (verbose) std::println("Finished MCMC.");
        delete(a);
        delete t->get_patindex();
        delete t->get_weights();
        delete(t);
        delete(p);
    }

    auto end = Clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::println("Time spent running: {} seconds",elapsed.count());
    if (verbose) std::println("End of ClonalFrame2 run.");
    return 0;
}