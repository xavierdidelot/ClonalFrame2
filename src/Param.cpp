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
    if (profiling) Profiler::prof.start(__func__);
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
    double sdRho=0.1;
    double sdDelta=10.0;
    double sdNu=1e-4;
    double acceptance_target=0.234;
    double acceptance_delta=1.0/acceptance_target/(1.0-acceptance_target);
    vector<int> accepted(4,0.0);
    std::normal_distribution<double> rnorm;

    for (int i=0;i<nbit;i++) {//Main MCMC loop
        if (verbose==false && floor((double)i/nbit*100)==(double)i/nbit*100) print_progress(i,nbit);
        if (verbose) vprint("Starting iteration {}",i);

        //Tree move
        Tree * tree2=Tree::copy_tree(tree);
        tree2->rSPR();
        tree2->jointml();
        double ll2=tree2->likelihood_anc_rec(rho,delta,nu);
        if (log(runif(rng))<ll2-ll) {
            accepted[0]++;
            if (verbose) vprint("Tree move accepted ll={:.2f} ll2={:.2f} paccept={:.2f}",ll,ll2,(double)accepted[0]/(i+1));
            ll=ll2;
            delete tree;
            tree=tree2;
        } else {
            if (verbose) vprint("Tree move rejected ll={:.2f} ll2={:.2f} paccept={:.2f}",ll,ll2,(double)accepted[0]/(i+1));
            delete tree2;
        }

        //#Adaptive MH move for rho using flat prior //Gamma(1e-3,1e3) prior
        rnorm.param(std::normal_distribution<double>::param_type(rho,sdRho));
        double rho2=abs(rnorm(rng));
        ll2=tree->likelihood_anc_rec(rho2,delta,nu);
        double mh=ll2-ll;//+dgamma(rho2,shape=1e-3,scale=1e3,log=T)-dgamma(rho,shape=1e-3,scale=1e3,log=T);
        if (log(runif(rng))<mh) {
            accepted[1]++;
            if (verbose) vprint("Rho move accepted rho={:.2f} rho2={:.2f} ll={:.2f} ll2={:.2f} sdRho={:.2f} paccept={:.2f}",rho,rho2,ll,ll2,sdRho,(double)accepted[1]/(i+1));
            ll=ll2;
            rho=rho2;
        } else {
            if (verbose) vprint("Rho move rejected rho={:.2f} rho2={:.2f} ll={:.2f} ll2={:.2f} sdRho={:.2f} paccept={:.2f}",rho,rho2,ll,ll2,sdRho,(double)accepted[1]/(i+1));
        }
        sdRho=sdRho*exp(acceptance_delta/(i+1)*(std::min(1.0,exp(mh))-acceptance_target));

        //#Adaptive MH move for delta using flat prior //Gamma(1e-3,1e3) prior
        rnorm.param(std::normal_distribution<double>::param_type(delta,sdDelta));
        double delta2=abs(rnorm(rng));
        ll2=tree->likelihood_anc_rec(rho,delta2,nu);
        mh=ll2-ll;//+dgamma(rho2,shape=1e-3,scale=1e3,log=T)-dgamma(rho,shape=1e-3,scale=1e3,log=T);
        if (log(runif(rng))<mh) {
            accepted[2]++;
            if (verbose) vprint("Delta move accepted delta={:.2f} delta2={:.2f} ll={:.2f} ll2={:.2f} sdDelta={:.2f} paccept={:.2f}",delta,delta2,ll,ll2,sdDelta,(double)accepted[2]/(i+1));
            ll=ll2;
            delta=delta2;
        } else {
            if (verbose) vprint("Delta move rejected delta={:.2f} delta2={:.2f} ll={:.2f} ll2={:.2f} sdDelta={:.2f} paccept={:.2f}",delta,delta2,ll,ll2,sdDelta,(double)accepted[2]/(i+1));
        }
        sdDelta=sdDelta*exp(acceptance_delta/(i+1)*(std::min(1.0,exp(mh))-acceptance_target));

        //#Adaptive MH move for nu using flat prior //Gamma(1e-3,1e3) prior
        rnorm.param(std::normal_distribution<double>::param_type(nu,sdNu));
        double nu2=abs(rnorm(rng));
        ll2=tree->likelihood_anc_rec(rho,delta,nu2);
        mh=ll2-ll;//+dgamma(rho2,shape=1e-3,scale=1e3,log=T)-dgamma(rho,shape=1e-3,scale=1e3,log=T);
        if (log(runif(rng))<mh) {
            accepted[3]++;
            if (verbose) vprint("Nu move accepted nu={:.2f} nu2={:.2f} ll={:.2f} ll2={:.2f} sdNu={:.2f} paccept={:.2f}",nu,nu2,ll,ll2,sdNu,(double)accepted[3]/(i+1));
            ll=ll2;
            nu=nu2;
        } else {
            if (verbose) vprint("Nu move rejected nu={:.2f} nu2={:.2f} ll={:.2f} ll2={:.2f} sdNu={:.2f} paccept={:.2f}",nu,nu2,ll,ll2,sdNu,(double)accepted[3]/(i+1));
        }
        sdNu=sdNu*exp(acceptance_delta/(i+1)*(std::min(1.0,exp(mh))-acceptance_target));

        //Writing to output file
        if (i>burn && i%thin==0) {
            std::println(file,"{}",tree->to_string());
        }
    }//Main MCMC loop
    if (verbose==false && nbit>0) print_progress(nbit,nbit);
    std::println(file,"");
    file.close();
    if (profiling) Profiler::prof.stop(__func__);
    if (verbose) vprint("Finished MCMC.");
}
