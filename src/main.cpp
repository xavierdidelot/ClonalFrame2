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
    -c NUM      Number of cores to use (default is all available)
    -S STR      Generate a simulated dataset and run on it
                STR is comma separated values for number of sequences,
                length of sequences, average branch length, rho, delta and nu
    -p          Activate profiler
    -v          Verbose mode
)";

int main(int argc, const char *argv[]) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto *orphan = new Node();orphan->set_dist_father(1.0); //just to test leak detection is working

    int burn = 10000, nbit = 20000, thin = 10;
    int i=1;
    bool simu=false;
    string simu_str;
    vector<string> remainder;
    while (i < argc) {
        string arg = argv[i];
        if (arg[0]=='-') {
            switch (arg[1]) {
                case('x'):burn =  std::stoi(argv[++i]) ;break;
                case('y'):nbit =  std::stoi(argv[++i]) ;break;
                case('z'):thin =  std::stoi(argv[++i]) ;break;
                case('s'):setSeed(std::stoi(argv[++i]));break;
                case('c'):n_cores =  std::stoi(argv[++i]);break;
                case('S'):simu = true;simu_str=argv[++i];break;
                case('v'):verbose = true;break;
                case('p'):profiling = true;break;
                default:std::println("Ignoring incorrect option -{}",arg[1]);
            }
        } else remainder.push_back(arg);
        i++;
    }

    vprint("This is ClonalFrame2, version {}", VERSION_STRING);
    if (remainder.size() != 2 && simu == false) {
        std::println("Need two non-option arguments but found {}", remainder.size());
        std::println("{}", help);
        return 0;
    }

    string cmd;
    for (i = 0; i < argc; i++) {cmd.append(argv[i]);cmd.append(" ");}
    std::println("Command used: {}", cmd);

    omp_set_num_threads(n_cores);
    std::println("Using {} cores.",n_cores);

    Alignment *a;
    string outfile;
    if (simu) {
        int n, l;
        double bralen, rho, delta, nu;
        std::sscanf(simu_str.c_str(), "%d,%d,%lf,%lf,%lf,%lf", &n, &l, &bralen, &rho, &delta, &nu);
        if (verbose) vprint("Generating true tree with {} leaves and mean branch length {}...", n, bralen);
        auto t0 = Tree::random_tree(n, bralen);
        if (verbose) vprint("Generating simulated alignment of length {} with rho={},delta={},nu={}...", l, rho, delta,
                            nu);
        t0->simulate_sequences_rec(l, rho, delta, nu);
        a = Alignment::make_from_tree(t0, false);
        delete t0->get_alignment();
        delete t0;
        outfile = "";
    } else {
        const auto inputfile = string(remainder[0]);
        outfile = string(remainder[1]);
        a = Alignment::read_fasta(inputfile);
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