//
// Created by Didelot, Xavier on 21/10/2025.
//

#include "Param.h"

Param::Param(Tree *t) {
    tree=t;
    rho=1;
    nu=0.01;
    delta=500.0;
    update_tree=1;
    update_rho=1;
    update_delta=1;
    update_nu=1;
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
    std::println(file, "{{");
    std::println(file, "\t\"parameters\": {{");
    std::println(file, "\t\t\"burn\": {},",burn);
    std::println(file, "\t\t\"nbit\": {},",nbit);
    std::println(file, "\t\t\"thin\": {}",thin);
    std::println(file, "\t}},");
    std::println(file, "\t\"iterations\": [");
    std::uniform_real_distribution runif(0.0,1.0);
    if (verbose) vprint("Performing initial ancestral state reconstruction...");
    tree->jointml();
    if (verbose) vprint("Performing initial likelihood computation...");
    ll=tree->likelihood_anc_rec(rho,delta,nu);
    tree_prior=tree->prior();
    double sdRho=0.1;
    double sdDelta=10.0;
    double sdNu=1e-4;
    double acceptance_target=0.234;
    double acceptance_delta=1.0/acceptance_target/(1.0-acceptance_target);
    double ll2;
    double mh;
    paccepted=vector(4,0.0);
    std::normal_distribution<double> rnorm;
    bool first_write=true;

    for (int i=0;i<nbit;i++) {//Main MCMC loop
        if (verbose==false && floor((double)i/nbit*100)==(double)i/nbit*100) print_progress(i,nbit);
        if (verbose) vprint("Starting iteration {}",i);

        //Writing to output file
        if (i%thin==0) {
            if (!first_write) {
                std::println(file, "\t,");
            } else first_write=false;
            to_json_file(&file);
        }

        //Tree move
        for (int rep=0;rep<update_tree;rep++) {
            Tree * tree2=Tree::copy_tree(tree);
            tree2->rSPR();
            tree2->jointml();
            ll2=tree2->likelihood_anc_rec(rho,delta,nu);
            double tree_prior2=tree2->prior();
            if (log(runif(rng))<ll2-ll+tree_prior2-tree_prior) {
                paccepted[0]=(paccepted[0]*i+1)/(i+1);
                if (verbose) vprint("Tree move accepted ll={:.2e} ll2={:.2e} paccept={:.2e}",ll,ll2,paccepted[0]);
                ll=ll2;
                tree_prior=tree_prior2;
                delete tree;
                tree=tree2;
            } else {
                paccepted[0]=(paccepted[0]*i)/(i+1);
                if (verbose) vprint("Tree move rejected ll={:.2e} ll2={:.2e} paccept={:.2e}",ll,ll2,paccepted[0]);
                delete tree2;
            }
        }

        //#Adaptive MH move for rho using flat prior //Gamma(1e-3,1e3) prior
        for (int rep=0;rep<update_rho;rep++) {
            rnorm.param(std::normal_distribution<double>::param_type(rho,sdRho));
            double rho2=std::fabs(rnorm(rng));
            ll2=tree->likelihood_anc_rec(rho2,delta,nu);
            mh=ll2-ll;//+dgamma(rho2,shape=1e-3,scale=1e3,log=T)-dgamma(rho,shape=1e-3,scale=1e3,log=T);
            if (log(runif(rng))<mh) {
                paccepted[1]=(paccepted[1]*i+1)/(i+1);
                if (verbose) vprint("Rho move accepted rho={:.2e} rho2={:.2e} ll={:.2e} ll2={:.2e} sdRho={:.2e} paccept={:.2f}",rho,rho2,ll,ll2,sdRho,paccepted[1]);
                ll=ll2;
                rho=rho2;
            } else {
                paccepted[1]=(paccepted[1]*i)/(i+1);
                if (verbose) vprint("Rho move rejected rho={:.2e} rho2={:.2e} ll={:.2e} ll2={:.2e} sdRho={:.2e} paccept={:.2f}",rho,rho2,ll,ll2,sdRho,paccepted[1]);
            }
            sdRho=sdRho*exp(acceptance_delta/(i+1)*(std::min(1.0,exp(mh))-acceptance_target));
        }

        //#Adaptive MH move for delta using flat prior //Gamma(1e-3,1e3) prior
        for (int rep=0;rep<update_delta;rep++) {
            rnorm.param(std::normal_distribution<double>::param_type(delta,sdDelta));
            double delta2=std::fabs(rnorm(rng));
            ll2=tree->likelihood_anc_rec(rho,delta2,nu);
            mh=ll2-ll;//+dgamma(rho2,shape=1e-3,scale=1e3,log=T)-dgamma(rho,shape=1e-3,scale=1e3,log=T);
            if (log(runif(rng))<mh) {
                paccepted[2]=(paccepted[2]*i+1)/(i+1);
                if (verbose) vprint("Delta move accepted delta={:.2e} delta2={:.2e} ll={:.2e} ll2={:.2e} sdDelta={:.2e} paccept={:.2f}",delta,delta2,ll,ll2,sdDelta,paccepted[2]);
                ll=ll2;
                delta=delta2;
            } else {
                paccepted[02]=(paccepted[2]*i)/(i+1);
                if (verbose) vprint("Delta move rejected delta={:.2e} delta2={:.2e} ll={:.2e} ll2={:.2e} sdDelta={:.2e} paccept={:.2f}",delta,delta2,ll,ll2,sdDelta,paccepted[2]);
            }
            sdDelta=sdDelta*exp(acceptance_delta/(i+1)*(std::min(1.0,exp(mh))-acceptance_target));
        }

        //#Adaptive MH move for nu using flat prior //Gamma(1e-3,1e3) prior
        for (int rep=0;rep<update_nu;rep++) {
            rnorm.param(std::normal_distribution<double>::param_type(nu,sdNu));
            double nu2=std::fabs(rnorm(rng));
            ll2=tree->likelihood_anc_rec(rho,delta,nu2);
            mh=ll2-ll;//+dgamma(rho2,shape=1e-3,scale=1e3,log=T)-dgamma(rho,shape=1e-3,scale=1e3,log=T);
            if (log(runif(rng))<mh) {
                paccepted[3]=(paccepted[3]*i+1)/(i+1);
                if (verbose) vprint("Nu move accepted nu={:.2e} nu2={:.2e} ll={:.2e} ll2={:.2e} sdNu={:.2e} paccept={:.2f}",nu,nu2,ll,ll2,sdNu,paccepted[3]);
                ll=ll2;
                nu=nu2;
            } else {
                paccepted[3]=(paccepted[3]*i)/(i+1);
                if (verbose) vprint("Nu move rejected nu={:.2e} nu2={:.2e} ll={:.2e} ll2={:.2e} sdNu={:.2e} paccept={:.2f}",nu,nu2,ll,ll2,sdNu,paccepted[3]);
            }
            sdNu=sdNu*exp(acceptance_delta/(i+1)*(std::min(1.0,exp(mh))-acceptance_target));
        }

    }//End of main MCMC loop
    if (verbose==false && nbit>0) print_progress(nbit,nbit);
    std::println(file, "\t]");
    std::println(file, "}}");
    file.close();
    if (profiling) Profiler::prof.stop(__func__);
    if (verbose) vprint("Finished MCMC.");
}

void Param::to_json_file(std::ofstream * file) {
    std::println(*file, "\t\t{{");
    std::println(*file,"\t\t\t\"likelihood\": {:.6e},",ll);
    std::println(*file,"\t\t\t\"tree_prior\": {:.6e},",tree_prior);
    std::println(*file,"\t\t\t\"sum_bra_len\": {:.6e},",tree->sum_bra_len());
    std::println(*file,"\t\t\t\"tree\": \"{}\",",tree->to_string());
    std::println(*file,"\t\t\t\"rho\": {:.6e},",rho);
    std::println(*file,"\t\t\t\"delta\": {:.6e},",delta);
    std::println(*file,"\t\t\t\"nu\": {:.6e},",nu);
    std::println(*file,"\t\t\t\"paccepted_tree\": {:.6e},",paccepted[0]);
    std::println(*file,"\t\t\t\"paccepted_rho\": {:.6e},",paccepted[1]);
    std::println(*file,"\t\t\t\"paccepted_delta\": {:.6e},",paccepted[2]);
    std::println(*file,"\t\t\t\"paccepted_nu\": {:.6e}",paccepted[3]);
    std::println(*file, "\t\t}}");
}

void Param::select_updates(const string &s) {
    std::sscanf(s.c_str(), "%d,%d,%d,%d", &update_tree, &update_rho, &update_delta, &update_nu);
}


