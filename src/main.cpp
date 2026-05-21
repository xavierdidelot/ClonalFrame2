#include <iostream>
#include <chrono>

#include "Alignment.h"
#include "Tree.h"
#include "Param.h"
#include "RNGGlobal.h"

constexpr const char* VERSION_STRING = "0.1";

static string helptxt =
        R"(Usage: ClonalFrame2 [OPTIONS]

    Options:
    -i FILE     Select the input file
    -o FILE     Select the output file
    -x NUM      Sets the number of burn-in iterations (default is 10000)
    -y NUM      Total number of iterations (default is 20000)
    -z NUM      Sets the number of iterations between samples (default is 10)
    -s NUM      Use given seed to initiate random number generator
    -c NUM      Number of cores to use (default is all available)
    -S STR      Generate a simulated dataset and run on it
                STR is comma separated values for number of sequences,
                length of sequences, average branch length, rho, delta and nu
    -p          Activate profiler
    -v          Verbose mode
    -h          Show help
)";

int main(int argc, const char *argv[]) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto *orphan = new Node();orphan->set_dist_father(1.0); //just to test leak detection is working

    int burn = 10000, nbit = 20000, thin = 10;
    int i=1;
    string infile="";
    string outfile="";
    string simu_str="";
    vector<string> remainder;
    bool help=false;
    while (i < argc) {
        string arg = argv[i];
        if (arg[0]=='-') {
            switch (arg[1]) {
                case('i'):infile=argv[++i];break;
                case('o'):outfile=argv[++i];break;
                case('x'):burn =  std::stoi(argv[++i]) ;break;
                case('y'):nbit =  std::stoi(argv[++i]) ;break;
                case('z'):thin =  std::stoi(argv[++i]) ;break;
                case('s'):setSeed(std::stoi(argv[++i]));break;
                case('c'):n_cores =  std::stoi(argv[++i]);break;
                case('S'):simu_str=argv[++i];break;
                case('v'):verbose = true;break;
                case('p'):profiling = true;break;
                case('h'):help=true;break;
                default:std::println("Ignoring incorrect option -{}",arg[1]);
            }
        } else remainder.push_back(arg);
        i++;
    }

    vprint("This is ClonalFrame2, version {}", VERSION_STRING);

    if (help) {
        std::println("{}", helptxt);
        return 0;
    }


    if (remainder.size() != 0) {
        std::println("Need zero non-option arguments but found {}", remainder.size());
        return 0;
    }

    if (simu_str.empty() && infile.empty()) {
        std::println("Need either -i or -S option.");
        return 0;
    }

    if (!simu_str.empty() && !infile.empty()) {
        std::println("Cannot have both -i and -S options.");
        return 0;
    }

    string cmd;
    for (i = 0; i < argc; i++) {cmd.append(argv[i]);cmd.append(" ");}
    std::println("Command used: {}", cmd);

    omp_set_num_threads(n_cores);
    std::println("Using {} cores.",n_cores);

    Alignment *a;
    if (simu_str.empty()) {
        a = Alignment::read_fasta(infile);
    } else {
        int n, l;
        double bralen, rho, delta, nu;
        std::sscanf(simu_str.c_str(), "%d,%d,%lf,%lf,%lf,%lf", &n, &l, &bralen, &rho, &delta, &nu);
        if (verbose) vprint("Generating true tree with {} leaves and mean branch length {}...", n, bralen);
        auto t0 = Tree::random_tree(n, bralen);
        if (verbose) vprint("Generating simulated alignment of length {} with rho={},delta={},nu={}...", l, rho, delta, nu);
        t0->simulate_sequences_rec(l, rho, delta, nu);
        a = t0->get_alignment();
        delete t0;
    }
    if (verbose) vprint("Building initial tree...");
    auto mat = a->distmat();
    auto t = Tree::nj(mat, a->get_names());
    delete mat;
    t->set_sequences(a);
    if (verbose) vprint("Performing initial ancestral reconstruction...");
    t->jointml();
    auto p = new Param(t);
    p->mcmc(burn, nbit, thin, outfile);
    t = p->get_tree();
    delete(a);
    delete(t);
    delete(p);

    auto end = Clock::now();
    std::chrono::duration<double> elapsed = end - start;
    vprint("Time spent running: {} seconds",elapsed.count());
    if (verbose) vprint("End of ClonalFrame2 run.");
    return 0;
}
