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
    int bar_width=60;
    float fraction = static_cast<float>(current) / total;
    int filled = static_cast<int>(fraction * bar_width);
    string bar(filled, '=');
    if (filled < bar_width) bar += '>';
    bar += string(bar_width - bar.size(), ' ');
    std::print("\r[{}] {}%  ", bar, static_cast<int>(fraction * 100));
    if (current == total) std::println("");;
}

void Param::mcmc(int burn, int nbit, int thin, string outfile) {
    std::ofstream file;
    if (!outfile.empty()) {
        file.open(outfile);
        if (file.is_open() == false) std::println(std::cerr,"Error opening file {}",outfile);
    }
    std::uniform_real_distribution runif(0.0,1.0);
    ll=tree->likelihood();

    for (int i=0;i<nbit;i++) {
        if (verbose==false && floor(i/nbit*100)==0) print_progress(i,nbit);
        if (verbose) std::println("Starting iteration {}",i);
        Tree * tree2=Tree::copy_tree(tree);
        tree2->rSPR();
        double ll2=tree2->likelihood();
        if (log(runif(rng))<ll2-ll) {
            //Accept move
            ll=ll2;
            delete tree;
            tree=tree2;
        } else {
            //Reject move
            delete tree2;
        }
        if (i>burn && i%thin==0) {
            std::println(file,"{}",tree->to_string());
        }
    }
    if (verbose==false) print_progress(nbit,nbit);
    std::println(file,"");
    file.close();
}
