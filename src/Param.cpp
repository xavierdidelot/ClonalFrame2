//
// Created by Didelot, Xavier on 21/10/2025.
//

#include "Param.h"

Param::Param(Tree *t) {
    tree=t;
    rho=1;
    nu=0.01;
    delta=500.0;
}

Param::~Param() {
}

void Param::print_progress(int current, int total) {
    int bar_width=50;
    float fraction = static_cast<float>(current) / total;
    int filled = static_cast<int>(fraction * bar_width);
    string bar(filled, '=');
    if (filled < bar_width) bar += '>';
    bar += string(bar_width - bar.size(), ' ');
    std::print("\r[{}] {}%  ", bar, static_cast<int>(fraction * 100));
    if (current == total) std::println("");;
    std::flush(std::cout);
}

void Param::mcmc(int burn, int nbit, int thin, string outfile) {
    if (verbose) vprint("Starting MCMC...");
    std::ofstream file;
    if (!outfile.empty()) {
        file.open(outfile);
        if (file.is_open() == false) std::println(std::cerr,"Error opening file {}",outfile);
    }
    std::uniform_real_distribution runif(0.0,1.0);
    if (verbose) vprint("Performing initial ancestral state reconstruction...");
    tree->jointml();
    if (verbose) vprint("Performing initial likelihood computation...");
    ll=tree->likelihood_anc_rec(rho,delta,nu);

    for (int i=0;i<nbit;i++) {
        if (verbose==false && floor((double)i/nbit*100)==(double)i/nbit*100) print_progress(i,nbit);
        if (verbose) vprint("Starting iteration {}",i);
        Tree * tree2=Tree::copy_tree(tree);
        tree2->rSPR();
        tree2->jointml();
        double ll2=tree2->likelihood_anc_rec(rho,delta,nu);
        if (log(runif(rng))<ll2-ll) {
            //Accept move
            if (verbose) vprint("Move accepted {} {}",ll,ll2);
            ll=ll2;
            delete tree;
            tree=tree2;
        } else {
            //Reject move
            if (verbose) vprint("Move rejected {} {}",ll,ll2);
            delete tree2;
        }
        if (i>burn && i%thin==0) {
            std::println(file,"{}",tree->to_string());
        }
    }
    if (verbose==false && nbit>0) print_progress(nbit,nbit);
    std::println(file,"");
    file.close();
    if (verbose) vprint("Finished MCMC.");
}
