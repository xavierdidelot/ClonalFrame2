//
// Created by Didelot, Xavier on 21/10/2025.
//

#include "Param.h"

Param::Param(Tree *t) {
    tree=t;
    rho=0.1;
    nu=0.05;
    delta=500.0;
    updates=vector<int>(5,1);
    prior_params={1,1,500,500,0.05,0.05,1e-4,1e-4};
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
    std::print(std::cout, "\r[{}] {}%  ", bar, static_cast<int>(fraction * 100));
    if (current == total) std::println(std::cout, "");
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
    prior_tree=tree->prior(this);
    prior_rho=dgamma_mean_sd(rho,prior_params[0],prior_params[1]);
    prior_delta=dgamma_mean_sd(delta,prior_params[2],prior_params[3]);
    prior_nu=dgamma_mean_sd(nu,prior_params[4],prior_params[5]);
    sd=vector<double> {0.1,10.0,1e-3,1e-3};
    adapt_nstart=vector(4,5.0*adapt_scale);
    adapt_lambda=vector(4,1.0);
    naccepted=vector(5,0);
    std::normal_distribution<double> rnorm;

    //Writing initial state to output file
    state_to_json_file(0);

    for (int i=0;i<nbit;i++) {//Main MCMC loop
        if (!verbose && nbit > 0 && i > 0 && (100 * i) / nbit != (100 * (i - 1)) / nbit) print_progress(i, nbit);
        if (verbose) vprint("Starting iteration {}",i);

        //#Adaptive MH move for rho using Gamma prior
        for (int rep=0;rep<updates[0];rep++) {
            rnorm.param(std::normal_distribution<double>::param_type(rho,sd[0]));
            double rho2=std::fabs(rnorm(rng));
            double prior_rho2=dgamma_mean_sd(rho2,prior_params[0],prior_params[1]);
            double ll2=tree->likelihood_anc_rec(rho2,delta,nu);
            double mh=ll2-ll+prior_rho2-prior_rho;
            if (log(runif(rng))<mh) {
                naccepted[0]++;
                if (verbose) vprint("Rho move accepted rho={:.2e} rho2={:.2e} ll={:.2e} ll2={:.2e} sd_rho={:.2e}",rho,rho2,ll,ll2,sd[0]);
                ll=ll2;
                rho=rho2;
                prior_rho=prior_rho2;
            } else {
                if (verbose) vprint("Rho move rejected rho={:.2e} rho2={:.2e} ll={:.2e} ll2={:.2e} sd_rho={:.2e}",rho,rho2,ll,ll2,sd[0]);
            }
            adapt(0,mh,i,rep);
        }

        //#Adaptive MH move for delta using Gamma prior
        for (int rep=0;rep<updates[1];rep++) {
            rnorm.param(std::normal_distribution<double>::param_type(delta,sd[1]));
            double delta2=std::fabs(rnorm(rng));
            double prior_delta2=dgamma_mean_sd(delta2,prior_params[2],prior_params[3]);
            double ll2=tree->likelihood_anc_rec(rho,delta2,nu);
            double mh=ll2-ll+prior_delta2-prior_delta;
            if (log(runif(rng))<mh) {
                naccepted[1]++;
                if (verbose) vprint("Delta move accepted delta={:.2e} delta2={:.2e} ll={:.2e} ll2={:.2e} sd_delta={:.2e}",delta,delta2,ll,ll2,sd[1]);
                ll=ll2;
                delta=delta2;
                prior_delta=prior_delta2;
            } else {
                if (verbose) vprint("Delta move rejected delta={:.2e} delta2={:.2e} ll={:.2e} ll2={:.2e} sd_delta={:.2e}",delta,delta2,ll,ll2,sd[1]);
            }
            adapt(1,mh,i,rep);
        }

        //#Adaptive MH move for nu using Gamma prior
        for (int rep=0;rep<updates[2];rep++) {
            rnorm.param(std::normal_distribution<double>::param_type(nu,sd[2]));
            double nu2=std::fabs(rnorm(rng));
            double prior_nu2=dgamma_mean_sd(nu2,prior_params[4],prior_params[5]);
            double ll2=tree->likelihood_anc_rec(rho,delta,nu2);
            double mh=ll2-ll+prior_nu2-prior_nu;
            if (log(runif(rng))<mh) {
                naccepted[2]++;
                if (verbose) vprint("Nu move accepted nu={:.2e} nu2={:.2e} ll={:.2e} ll2={:.2e} sd_nu={:.2e}",nu,nu2,ll,ll2,sd[2]);
                ll=ll2;
                nu=nu2;
                prior_nu=prior_nu2;
            } else {
                if (verbose) vprint("Nu move rejected nu={:.2e} nu2={:.2e} ll={:.2e} ll2={:.2e} sd_nu={:.2e}",nu,nu2,ll,ll2,sd[2]);
            }
            adapt(2,mh,i,rep);
        }

        //Adaptive MH move for branch lengths
        for (int rep=0;rep<updates[3];rep++) {
            Tree * tree2=Tree::copy_tree(tree);
            rnorm.param(std::normal_distribution<double>::param_type(0,sd[3]));
            double lqr=0.0;
            for (size_t b=0;b<2*tree2->get_n()-3;b++) {
                double cur=tree2->get_node(b)->get_dist_father();
                double pro=cur*std::exp(rnorm(rng));
                tree2->get_node(b)->set_dist_father(pro);
                lqr=lqr+std::log(pro/cur);
            }
            tree2->jointml();
            double ll2=tree2->likelihood_anc_rec(rho,delta,nu);
            double prior_tree2=tree2->prior(this);
            double mh=ll2-ll+prior_tree2-prior_tree+lqr;
            if (log(runif(rng))<mh) {
                naccepted[3]++;
                if (verbose) vprint("Branch lengths move accepted ll={:.2e} ll2={:.2e} lp={:.2e} lp2={:.2e} mh={:2e}",ll,ll2,prior_tree,prior_tree2,mh);
                ll=ll2;
                prior_tree=prior_tree2;
                delete tree;
                tree=tree2;
            } else {
                if (verbose) vprint("Branch lengths move rejected ll={:.2e} ll2={:.2e} lp={:.2e} lp2={:.2e} mh={:2e}",ll,ll2,prior_tree,prior_tree2,mh);
                delete tree2;
            }
            adapt(3,mh,i,rep);
        }

        //SPR move
        for (int rep=0;rep<updates[4];rep++) {
            Tree * tree2=Tree::copy_tree(tree);
            string name_moving;
            double q12=tree2->rSPR(this,&name_moving);
            tree2->jointml();
            double q21=tree2->lq_rSPR(tree,this,name_moving);
            double ll2=tree2->likelihood_anc_rec(rho,delta,nu);
            double prior_tree2=tree2->prior(this);
            double mh=ll2-ll+prior_tree2-prior_tree+q21-q12;
            if (log(runif(rng))<mh) {
                naccepted[4]++;
                if (verbose) vprint("Tree move accepted ll2-ll={:.2e} lp2-lp={:.2e} lq21-lq12={:.2e} mh={:.2e}",ll2-ll,prior_tree2-prior_tree,q21-q12,mh);
                ll=ll2;
                prior_tree=prior_tree2;
                delete tree;
                tree=tree2;
            } else {
                if (verbose) vprint("Tree move rejected ll2-ll={:.2e} lp2-lp={:.2e} lq21-lq12={:.2e} mh={:.2e}",ll2-ll,prior_tree2-prior_tree,q21-q12,mh);
                delete tree2;
            }
        }

        //Writing to output file
        if ((i+1)%thin==0) {
            std::println(out, "\t,");
            state_to_json_file(i+1);
        }

    }//End of main MCMC loop
    if (verbose==false && nbit>0) print_progress(nbit,nbit);
    std::println(out, "\t]");
    if (profiling) Profiler::prof.stop(__func__);
    if (verbose) vprint("Finished MCMC.");
}

void Param::adapt(int index,double mh,int i,int rep) {
    double prob_accept=mh>=0?1.0:std::exp(mh);//this is same as min(1,exp(mh))
    int it=1+rep+i*updates[index];
    sd[index]*=exp(0.5*adapt_scale/(it+adapt_nstart[index])*(prob_accept-adapt_target));
    double lambda=sd[index]*sd[index];
    if (std::abs(std::log(lambda)-std::log(adapt_lambda[index]))>std::log(3.0)) {
        //Restart Robbins-Munro
        adapt_lambda[index]=lambda;
        adapt_nstart[index]=5.0*adapt_scale-it;
    }
}

void Param::state_to_json_file(int iter) {
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
    std::println(out,"\t\t\t\"sd_rho\": {:.6e},",sd[0]);
    std::println(out,"\t\t\t\"sd_delta\": {:.6e},",sd[1]);
    std::println(out,"\t\t\t\"sd_nu\": {:.6e},",sd[2]);
    std::println(out,"\t\t\t\"sd_lengths\": {:.6e},",sd[3]);
    std::println(out,"\t\t\t\"paccepted_rho\": {:.6e},",    iter*updates[0]==0?0.0:1.0*naccepted[0]/(iter*updates[0]));
    std::println(out,"\t\t\t\"paccepted_delta\": {:.6e},",  iter*updates[1]==0?0.0:1.0*naccepted[1]/(iter*updates[1]));
    std::println(out,"\t\t\t\"paccepted_nu\": {:.6e},",     iter*updates[2]==0?0.0:1.0*naccepted[2]/(iter*updates[2]));
    std::println(out,"\t\t\t\"paccepted_lengths\": {:.6e},",iter*updates[3]==0?0.0:1.0*naccepted[3]/(iter*updates[3]));
    std::println(out,"\t\t\t\"paccepted_spr\": {:.6e}",     iter*updates[4]==0?0.0:1.0*naccepted[4]/(iter*updates[4]));
    std::println(out, "\t\t}}");
}
