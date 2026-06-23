#include <iostream>
#include <chrono>

#include "Alignment.h"
#include "Tree.h"
#include "Param.h"
#include "global.h"

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
    -b NUM      Mean of the prior on branch length
    -S STR      Generate a simulated dataset and run on it
                STR is 7 comma-separated values for number of sequences,
                length of sequences, mean branch length, rho, delta, nu
                and whether to start from the correct parameters
    -U STR      Indicates number of updates to use per iteration for each move
                STR is 5 int values for rho, delta, nu, lengths and SPR
    -p          Activate profiler
    -v          Verbose mode
    -h          Show help
)";

int main(int argc, const char *argv[]) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto *orphan = new Node();orphan->set_dist_father(1.0); //just to test leak detection is working

    int burn = 10000, nbit = 20000, thin = 10;
    bra_len_prior_mean=1e-3;
    int i=1;
    string infile="";
    string outfile="";
    string simu_str="";
    string update_str="1,1,1,1,1";
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
                case('b'):bra_len_prior_mean=std::stof(argv[++i]);break;
                case('S'):simu_str=argv[++i];break;
                case('U'):update_str=argv[++i];break;
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

    if (!outfile.empty()) {
        out.open(outfile);
        if (out.is_open() == false) std::println(std::cerr,"Unable to open output file {}",outfile);
    }
    std::println(out, "{{");

    string cmd;
    for (i = 0; i < argc; i++) {cmd.append(argv[i]);cmd.append(" ");}
    std::println("Command used: {}", cmd);
    std::println(out,"\t\"Command\": \"{}\",",cmd);

    omp_set_num_threads(n_cores);
    std::println("Using {} cores.",n_cores);

    Alignment *a;
    Tree *t=nullptr;
    Param *p;
    double bralen, rho, delta, nu;
    int start_correct=0;
    if (simu_str.empty()) {
        a = Alignment::read_fasta(infile);
        if (a->get_n()<3) {
            std::println("Need a minimum of three genomes.");
            return 0;
        }
    } else {
        int n, l;
        std::sscanf(simu_str.c_str(), "%d,%d,%lf,%lf,%lf,%lf,%d", &n, &l, &bralen, &rho, &delta, &nu, &start_correct);
        if (n<3) {
            std::println("Need a minimum of three genomes.");
            return 0;
        }
        if (verbose) vprint("Generating true tree with {} leaves and mean branch length {}...", n, bralen);
        t = Tree::random_tree(n, bralen);
        if (verbose) vprint("Generating simulated alignment of length {} with rho={},delta={},nu={}...", l, rho, delta, nu);
        std::println(out,"\t\"simulation_parameters\": {{");
        std::println(out,"\t\t\"tree\": \"{}\",",t->to_string());
        std::println(out,"\t\t\"length\": {},",l);
        std::println(out,"\t\t\"rho\": {:.6e},",rho);
        std::println(out,"\t\t\"delta\": {:.6e},",delta);
        std::println(out,"\t\t\"nu\": {:.6e},",nu);
        t->simulate_sequences_rec(l, rho, delta, nu);
        std::println(out,"\t}},");
        a = t->get_alignment();
    }

    if (start_correct==1) {
        t->set_sequences(a);
        if (verbose) vprint("Starting from correct tree and parameters...");
        p=new Param(t);
        p->set_rho(rho);
        p->set_delta(delta);
        p->set_nu(nu);
    } else {
        delete t;
        if (verbose) vprint("Building initial tree...");
        auto mat = a->distmat();
        t = Tree::nj(mat, a->get_names());
        delete mat;
        t->set_sequences(a);
        p = new Param(t);
    }
    if (verbose) vprint("Performing initial ancestral reconstruction...");
    t->jointml();
    p->select_updates(update_str);
    p->mcmc(burn, nbit, thin);
    t = p->get_tree();
    delete(a);
    delete(t);
    delete(p);
    std::println(out, "}}");
    out.close();

    auto end = Clock::now();
    std::chrono::duration<double> elapsed = end - start;
    vprint("Time spent running: {} seconds",elapsed.count());
    if (verbose) vprint("End of ClonalFrame2 run.");
    return 0;
}
