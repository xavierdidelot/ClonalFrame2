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
    -x NUM      Sets the total number of iterations (default is 10000)
    -y NUM      Sets the number of iterations between samples (default is x/2000)
    -s NUM      Use given seed to initiate random number generator
    -c NUM      Number of cores to use (default is all available)
    -S STR      Generate a simulated dataset and run on it
                STR is 7 comma-separated values for number of sequences,
                length of sequences, rho, delta, nu, branch length mean and sd
    -I STR      Specify starting parameter values (default is 1,500,0.05)
                STR is 3 comma-separated values for rho, delta and nu
                If -S is used without -I then start from correct values and tree
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
    #ifdef TESTLEAK
        //to test leak detection is working, need to compile with -DTESTLEAK for this to be included
        auto *orphan = new Node();orphan->set_dist_father(1.0);
    #endif
    int nbit = 20000, thin = -1;
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
                case('x'):nbit =  std::stoi(argv[++i]) ;break;
                case('y'):thin =  std::stoi(argv[++i]) ;break;
                case('s'):setSeed(std::stoi(argv[++i]));break;
                case('c'):n_cores =  std::stoi(argv[++i]);break;
                case('S'):simu_str=argv[++i];break;
                case('I'):starting_str=argv[++i];break;
                case('U'):update_str=argv[++i];break;
                case('P'):prior_params_str=argv[++i];break;
                case('v'):verbose = true;break;
                case('p'):profiling = true;break;
                case('h'):help=true;break;
                default:std::cout<<"Ignoring incorrect option -"<<arg[1]<<std::endl;
            }
        } else remainder.push_back(arg);
        i++;
    }
    if (thin==-1) thin=nbit/2000;
    if (thin==0) thin=1;

    std::cout<<vtime()<<"This is ClonalFrame2, version "<<VERSION_STRING<<std::endl;

    if (help) {
        std::cout<<helptxt<<std::endl;
        return 0;
    }


    if (remainder.size() != 0) {
        std::cerr<<"Need zero non-option arguments but found "<< remainder.size()<<std::endl;
        return 0;
    }

    if (simu_str.empty() && infile.empty()) {
        std::cerr<<"Need either -i or -S option."<<std::endl;
        return 0;
    }

    if (!simu_str.empty() && !infile.empty()) {
        std::cerr<<"Cannot have both -i and -S options."<<std::endl;
        return 0;
    }

    if (!outfile.empty()) {
        out.open(outfile);
        if (out.is_open() == false) std::cerr<<"Unable to open output file "<<outfile<<std::endl;
    }
    out<<"{\n";
    string cmd;
    for (i = 0; i < argc; i++) {cmd.append(argv[i]);cmd.append(" ");}
    std::cout<<vtime()<<"Command used: "<< cmd<<std::endl;
    out<<"\t\"command\": \""<<cmd<<"\",\n";

    n_cores=set_num_threads(n_cores);
    std::cout<<vtime()<<"Using "<<n_cores<<" cores."<<std::endl;

    Alignment *a;
    Tree *t=nullptr;
    Param *p;
    double b_mean, b_sd, rho, delta, nu;
    int start_correct=0;
    if (simu_str.empty()) {
        if (starting_str.empty()) starting_str="0.1,500,0.05";
        a = Alignment::read_fasta(infile);
        if (a->get_n()<3) {
            std::cerr<<"Need a minimum of three genomes."<<std::endl;
            return 0;
        }
    } else {
        int n, l;
        std::sscanf(simu_str.c_str(), "%d,%d,%lf,%lf,%lf,%lf,%lf", &n, &l, &rho, &delta, &nu, &b_mean, &b_sd);
        if (starting_str.empty()) start_correct=1;
        if (n<3) {
            std::cerr<<"Need a minimum of three genomes."<<std::endl;
            return 0;
        }
        if (verbose) std::cout<<vtime()<<"Generating true tree with "<<n<<" leaves and branch length with mean "<<b_mean<<" and sd "<<b_sd<<"..."<<std::endl;
        t = Tree::random_tree(n, b_mean, b_sd);
        if (verbose) std::cout<<vtime()<<"Generating simulated alignment of length "<<l<<" with rho="<<rho<<",delta="<<delta<<",nu="<<nu<<"..."<<std::endl;
        out << "\t\"simulation_parameters\": {\n"
            << "\t\t\"tree\": \"" << t->to_string() << "\",\n"
            << "\t\t\"length\": " << l << ",\n"
            << std::scientific << std::setprecision(6)
            << "\t\t\"rho\": " << rho << ",\n"
            << "\t\t\"delta\": " << delta << ",\n"
            << "\t\t\"nu\": " << nu << ",\n";
        t->simulate_sequences_rec(l, rho, delta, nu);
        out << "\t},\n";
        a = t->get_alignment();
    }

    if (start_correct==1) {
        t->set_sequences(a);
        if (verbose) std::cout<<vtime()<<"Starting from correct tree and parameters..."<<std::endl;
        p=new Param(t);
        if (rho>0) p->set_rho(rho); else p->set_rho(1e-5);
        if (delta>0) p->set_delta(delta); else p->set_delta(1e-5);
        if (nu>0) p->set_nu(nu); else p->set_nu(1e-5);
    } else {
        delete t;
        if (verbose) std::cout<<vtime()<<"Building initial tree..."<<std::endl;
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
    vector<int> * updates=p->get_updates();
    std::sscanf(update_str.c_str(), "%d,%d,%d,%d,%d", &updates->at(0),&updates->at(1),&updates->at(2),&updates->at(3),&updates->at(4));
    p->mcmc(nbit, thin);
    t = p->get_tree();
    delete(a);
    delete(t);
    delete(p);
    out << "}\n";
    out.close();

    auto end = Clock::now();
    std::chrono::duration<double> elapsed = end - start;
    double time=elapsed.count();
    string unit="seconds";
    if (time>60) {time=time/60;unit="minutes";}
    if (time>60 && unit=="minutes") {time=time/60;unit="hours";}
    if (time>24 && unit=="hours") {time=time/24;unit="days";}
    std::cout<<vtime()<<"Time spent running: "<<std::setprecision(4)<<time<<" "<<unit<<std::endl;
    if (verbose) std::cout<<vtime()<<"End of ClonalFrame2 run."<<std::endl;
    return 0;
}
