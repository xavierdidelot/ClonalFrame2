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
    -x NUM      Sets the number of burn-in iterations (default is y/2)
    -y NUM      Total number of iterations (default is 10000)
    -z NUM      Sets the number of iterations between samples (default is y/2000)
    -s NUM      Use given seed to initiate random number generator
    -c NUM      Number of cores to use (default is all available)
    -S STR      Generate a simulated dataset and run on it
                STR is 7 comma-separated values for number of sequences,
                length of sequences, rho, delta, nu, branch length mean and sd
    -I STR      Specify starting values (default is 1,500,0.05,nj
                unless -S was used in which case default is correct values)
                STR is 4 comma-separated values for rho, delta, nu and tree
    -U STR      Indicates number of updates to use per iteration for each move
                STR is 5 int values for rho, delta, nu, lengths and SPR
    -P STR      Indicates prior as 8 comma-separated values for mean and sd
                of rho, delta, nu and b (default 1,1,500,500,0.05,0.05,1e-4,1e-4)
    -p          Activate profiler
    -v          Verbose mode
    -h          Show help
)";

int main(int argc, const char *argv[]) {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto *orphan = new Node();orphan->set_dist_father(1.0); //just to test leak detection is working

    int burn = -1, nbit = 20000, thin = -1;
    int i=1;
    string infile="";
    string outfile="";
    string simu_str="";
    string starting_str="";
    string prior_params_str="";
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
                case('S'):simu_str=argv[++i];break;
                case('I'):starting_str=argv[++i];break;
                case('U'):update_str=argv[++i];break;
                case('P'):prior_params_str=argv[++i];break;
                case('v'):verbose = true;break;
                case('p'):profiling = true;break;
                case('h'):help=true;break;
                default:std::println("Ignoring incorrect option -{}",arg[1]);
            }
        } else remainder.push_back(arg);
        i++;
    }
    if (burn==-1) burn=nbit/2;
    if (thin==-1) thin=nbit/2000;
    if (thin==0) thin=1;

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
    std::println(out,"\t\"command\": \"{}\",",cmd);

    omp_set_num_threads(n_cores);
    std::println("Using {} cores.",n_cores);

    Alignment *a;
    Tree *t=nullptr;
    Param *p;
    double b_mean, b_sd, rho, delta, nu;
    int start_correct=0;
    if (simu_str.empty()) {
        if (starting_str.empty()) starting_str="1,500,0.05,nj";
        a = Alignment::read_fasta(infile);
        if (a->get_n()<3) {
            std::println("Need a minimum of three genomes.");
            return 0;
        }
    } else {
        int n, l;
        std::sscanf(simu_str.c_str(), "%d,%d,%lf,%lf,%lf,%lf,%lf", &n, &l, &rho, &delta, &nu, &b_mean, &b_sd);
        if (starting_str.empty()) start_correct=1;
        if (n<3) {
            std::println("Need a minimum of three genomes.");
            return 0;
        }
        if (verbose) vprint("Generating true tree with {} leaves and branch length with mean {} and sd {}...", n, b_mean, b_sd);
        t = Tree::random_tree(n, b_mean, b_sd);
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
        if (rho>0) p->set_rho(rho); else p->set_rho(1e-5);
        if (delta>0) p->set_delta(delta); else p->set_delta(1e-5);
        if (nu>0) p->set_nu(nu); else p->set_nu(1e-5);
    } else {
        delete t;
        if (verbose) vprint("Building initial tree...");
        auto mat = a->distmat();
        t = Tree::nj(mat, a->get_names());
        delete mat;
        t->set_sequences(a);
        p = new Param(t);
        std::sscanf(starting_str.c_str(), "%lf,%lf,%lf", &rho, &delta, &nu);
        p->set_rho(rho);
        p->set_delta(delta);
        p->set_nu(nu);
    }
    if (prior_params_str.length()>0) {
        vector<double> * prior_params=p->get_prior_params();
        std::sscanf(prior_params_str.c_str(), "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &prior_params->at(0), &prior_params->at(1),
            &prior_params->at(2),&prior_params->at(3),&prior_params->at(4),&prior_params->at(5),&prior_params->at(6),&prior_params->at(7));
    }
    if (verbose) vprint("Performing initial ancestral reconstruction...");
    t->jointml();
    vector<int> * updates=p->get_updates();
    std::sscanf(update_str.c_str(), "%d,%d,%d,%d,%d", &updates->at(0),&updates->at(1),&updates->at(2),&updates->at(3),&updates->at(4));
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
