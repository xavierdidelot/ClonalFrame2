//
// Created by Didelot, Xavier on 21/10/2025.
//

#include "Param.h"

Param::Param(Tree *t) {
    tree=t;
    rho=1;
    nu=0.01;
    delta=500.0;
    update_rho=1;
    update_delta=1;
    update_nu=1;
    update_lengths=1;
    update_spr=1;
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

void Param::mcmc(int burn, int nbit, int thin) {
    if (verbose) vprint("Starting MCMC...");
    if (profiling) Profiler::prof.start(__func__);
    std::println(out, "\t\"mcmc_parameters\": {{");
    std::println(out, "\t\t\"burn\": {},",burn);
    std::println(out, "\t\t\"nbit\": {},",nbit);
    std::println(out, "\t\t\"thin\": {}",thin);
    std::println(out, "\t}},");
    std::uniform_real_distribution runif(0.0,1.0);
    if (verbose) vprint("Performing initial ancestral state reconstruction...");
    tree->jointml();
    tree->ancestral_to_json_file();
    std::println(out, "\t\"iterations\": [");
    if (verbose) vprint("Performing initial likelihood computation...");
    ll=tree->likelihood_anc_rec(rho,delta,nu);
    prior_tree=tree->prior();
    prior_rho=dgamma_mean_sd(rho,prior_params[0],prior_params[1]);
    prior_delta=dgamma_mean_sd(delta,prior_params[2],prior_params[3]);
    prior_nu=dgamma_mean_sd(nu,prior_params[4],prior_params[5]);
    double sdRho=0.1;
    double sdDelta=10.0;
    double sdNu=1e-3;
    double sdLengths=prior_params[6]*0.01;
    double acceptance_target=0.234;
    double acceptance_delta=1.0/acceptance_target/(1.0-acceptance_target);
    double ll2;
    double mh;
    paccepted=vector(5,0.0);
    std::normal_distribution<double> rnorm;
    bool first_write=true;

    for (int i=0;i<nbit;i++) {//Main MCMC loop
        if (verbose==false && floor((double)i/nbit*100)==(double)i/nbit*100) print_progress(i,nbit);
        if (verbose) vprint("Starting iteration {}",i);

        //Writing to output file
        if (i%thin==0) {
            if (!first_write) {
                std::println(out, "\t,");
            } else first_write=false;
            state_to_json_file();
        }

        //#Adaptive MH move for rho using Gamma prior
        for (int rep=0;rep<update_rho;rep++) {
            rnorm.param(std::normal_distribution<double>::param_type(rho,sdRho));
            double rho2=std::fabs(rnorm(rng));
            double prior_rho2=dgamma_mean_sd(rho2,prior_params[0],prior_params[1]);
            ll2=tree->likelihood_anc_rec(rho2,delta,nu);
            mh=ll2-ll+prior_rho2-prior_rho;
            if (log(runif(rng))<mh) {
                paccepted[0]=(paccepted[0]*(i*update_rho+rep)+1)/(i*update_rho+rep+1);
                if (verbose) vprint("Rho move accepted rho={:.2e} rho2={:.2e} ll={:.2e} ll2={:.2e} sdRho={:.2e} paccept={:.2f}",rho,rho2,ll,ll2,sdRho,paccepted[0]);
                ll=ll2;
                rho=rho2;
                prior_rho=prior_rho2;
            } else {
                paccepted[0]=(paccepted[0]*(i*update_rho+rep)+0)/(i*update_rho+rep+1);
                if (verbose) vprint("Rho move rejected rho={:.2e} rho2={:.2e} ll={:.2e} ll2={:.2e} sdRho={:.2e} paccept={:.2f}",rho,rho2,ll,ll2,sdRho,paccepted[0]);
            }
            sdRho=sdRho*exp(acceptance_delta/(i+1)*(std::min(1.0,exp(mh))-acceptance_target));
        }

        //#Adaptive MH move for delta using Gamma prior
        for (int rep=0;rep<update_delta;rep++) {
            rnorm.param(std::normal_distribution<double>::param_type(delta,sdDelta));
            double delta2=std::fabs(rnorm(rng));
            double prior_delta2=dgamma_mean_sd(delta2,prior_params[2],prior_params[3]);
            ll2=tree->likelihood_anc_rec(rho,delta2,nu);
            mh=ll2-ll+prior_delta2-prior_delta;
            if (log(runif(rng))<mh) {
                paccepted[1]=(paccepted[1]*(i*update_delta+rep)+1)/(i*update_delta+rep+1);
                if (verbose) vprint("Delta move accepted delta={:.2e} delta2={:.2e} ll={:.2e} ll2={:.2e} sdDelta={:.2e} paccept={:.2f}",delta,delta2,ll,ll2,sdDelta,paccepted[1]);
                ll=ll2;
                delta=delta2;
                prior_delta=prior_delta2;
            } else {
                paccepted[1]=(paccepted[1]*(i*update_delta+rep)+0)/(i*update_delta+rep+1);
                if (verbose) vprint("Delta move rejected delta={:.2e} delta2={:.2e} ll={:.2e} ll2={:.2e} sdDelta={:.2e} paccept={:.2f}",delta,delta2,ll,ll2,sdDelta,paccepted[1]);
            }
            sdDelta=sdDelta*exp(acceptance_delta/(i+1)*(std::min(1.0,exp(mh))-acceptance_target));
        }

        //#Adaptive MH move for nu using using Gamma prior
        for (int rep=0;rep<update_nu;rep++) {
            rnorm.param(std::normal_distribution<double>::param_type(nu,sdNu));
            double nu2=std::fabs(rnorm(rng));
            double prior_nu2=dgamma_mean_sd(nu2,prior_params[4],prior_params[5]);
            ll2=tree->likelihood_anc_rec(rho,delta,nu2);
            mh=ll2-ll+prior_nu2-prior_nu;
            if (log(runif(rng))<mh) {
                paccepted[2]=(paccepted[2]*(i*update_nu+rep)+1)/(i*update_nu+rep+1);
                if (verbose) vprint("Nu move accepted nu={:.2e} nu2={:.2e} ll={:.2e} ll2={:.2e} sdNu={:.2e} paccept={:.2f}",nu,nu2,ll,ll2,sdNu,paccepted[2]);
                ll=ll2;
                nu=nu2;
                prior_nu=prior_nu2;
            } else {
                paccepted[2]=(paccepted[2]*(i*update_nu+rep)+0)/(i*update_nu+rep+1);
                if (verbose) vprint("Nu move rejected nu={:.2e} nu2={:.2e} ll={:.2e} ll2={:.2e} sdNu={:.2e} paccept={:.2f}",nu,nu2,ll,ll2,sdNu,paccepted[2]);
            }
            sdNu=sdNu*exp(acceptance_delta/(i+1)*(std::min(1.0,exp(mh))-acceptance_target));
        }

        //Adaptive MH move for branch lengths
        for (int rep=0;rep<update_lengths;rep++) {
            Tree * tree2=Tree::copy_tree(tree);
            rnorm.param(std::normal_distribution<double>::param_type(0,sdLengths));
            for (size_t b=0;b<2*tree2->get_n()-3;b++) {
                tree2->get_node(b)->set_dist_father(std::fabs(tree2->get_node(b)->get_dist_father()+rnorm(rng)));
            }
            tree2->jointml();
            ll2=tree2->likelihood_anc_rec(rho,delta,nu);
            double prior_tree2=tree2->prior();
            mh=ll2-ll+prior_tree2-prior_tree;
            if (log(runif(rng))<mh) {
                paccepted[3]=(paccepted[3]*(i*update_lengths+rep)+1)/(i*update_lengths+rep+1);
                if (verbose) vprint("Branch lengths move accepted ll={:.2e} ll2={:.2e} lp={:.2e} lp2={:.2e} paccept={:.2e}",ll,ll2,prior_tree,prior_tree2,paccepted[3]);
                ll=ll2;
                prior_tree=prior_tree2;
                delete tree;
                tree=tree2;
            } else {
                paccepted[3]=(paccepted[3]*(i*update_lengths+rep)+0)/(i*update_lengths+rep+1);
                if (verbose) vprint("Branch lengths move rejected ll={:.2e} ll2={:.2e} lp={:.2e} lp2={:.2e} paccept={:.2e}",ll,ll2,prior_tree,prior_tree2,paccepted[3]);
                delete tree2;
            }
            sdLengths=sdLengths*exp(acceptance_delta/(i+1)*(std::min(1.0,exp(mh))-acceptance_target));

        }

        //SPR move
        for (int rep=0;rep<update_spr;rep++) {
            Tree * tree2=Tree::copy_tree(tree);
            string name_moving;
            double q12=tree2->rSPR(&name_moving);
            tree2->jointml();
            double q21=tree2->lq_rSPR(tree,name_moving);
            ll2=tree2->likelihood_anc_rec(rho,delta,nu);
            double prior_tree2=tree2->prior();
            if (log(runif(rng))<ll2-ll+prior_tree2-prior_tree+q21-q12) {
                paccepted[4]=(paccepted[4]*(i*update_spr+rep)+1)/(i*update_spr+rep+1);
                if (verbose) vprint("Tree move accepted ll={:.2e} ll2={:.2e} lp={:.2e} lp2={:.2e} lq21={:.2e} lq21={:.2e} paccept={:.2e}",ll,ll2,prior_tree,prior_tree2,q21,q12,paccepted[4]);
                ll=ll2;
                prior_tree=prior_tree2;
                delete tree;
                tree=tree2;
            } else {
                paccepted[4]=(paccepted[4]*(i*update_spr+rep)+0)/(i*update_spr+rep+1);
                if (verbose) vprint("Tree move rejected ll={:.2e} ll2={:.2e} lp={:.2e} lp2={:.2e} lq21={:.2e} lq21={:.2e} paccept={:.2e}",ll,ll2,prior_tree,prior_tree2,q21,q12,paccepted[4]);
                delete tree2;
            }
        }

    }//End of main MCMC loop
    if (verbose==false && nbit>0) print_progress(nbit,nbit);
    std::println(out, "\t]");
    if (profiling) Profiler::prof.stop(__func__);
    if (verbose) vprint("Finished MCMC.");
}

void Param::state_to_json_file() {
    std::println(out, "\t\t{{");
    std::println(out,"\t\t\t\"likelihood\": {:.6e},",ll);
    std::println(out,"\t\t\t\"prior_tree\": {:.6e},",prior_tree);
    std::println(out,"\t\t\t\"prior_rho\": {:.6e},",prior_rho);
    std::println(out,"\t\t\t\"prior_delta\": {:.6e},",prior_delta);
    std::println(out,"\t\t\t\"prior_nu\": {:.6e},",prior_nu);
    std::println(out,"\t\t\t\"sum_bra_len\": {:.6e},",tree->sum_bra_len());
    std::println(out,"\t\t\t\"tree\": \"{}\",",tree->to_string());
    std::println(out,"\t\t\t\"rho\": {:.6e},",rho);
    std::println(out,"\t\t\t\"delta\": {:.6e},",delta);
    std::println(out,"\t\t\t\"nu\": {:.6e},",nu);
    std::println(out,"\t\t\t\"paccepted_rho\": {:.6e},",paccepted[0]);
    std::println(out,"\t\t\t\"paccepted_delta\": {:.6e},",paccepted[1]);
    std::println(out,"\t\t\t\"paccepted_nu\": {:.6e},",paccepted[2]);
    std::println(out,"\t\t\t\"paccepted_lengths\": {:.6e},",paccepted[3]);
    std::println(out,"\t\t\t\"paccepted_spr\": {:.6e}",paccepted[4]);
    std::println(out, "\t\t}}");
}

void Param::select_updates(const string &s) {
    std::sscanf(s.c_str(), "%d,%d,%d,%d,%d", &update_rho, &update_delta, &update_nu,&update_lengths,&update_spr);
}

