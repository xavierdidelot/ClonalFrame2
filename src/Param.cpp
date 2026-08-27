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
    std::cout << '\r' << '[' << bar << "] " << static_cast<int>(fraction * 100) << "%  ";
    if (current == total) std::cout << "" << std::endl;
    std::flush(std::cout);
}

void Param::mcmc(int nbit, int thin) {
    if (verbose) std::cout<<vtime()<<"Starting MCMC..."<<std::endl;
    if (profiling) Profiler::prof.start(__func__);
    out << "\t\"mcmc_parameters\": {\n";
    out << "\t\t\"nbit\": "<<nbit<<",\n";
    out << "\t\t\"thin\": "<<thin<<"\n";
    out << "\t},\n";
    std::uniform_real_distribution runif(0.0,1.0);
    if (verbose) std::cout<<vtime()<<"Performing initial ancestral state reconstruction..."<<std::endl;
    tree->jointml();
    tree->ancestral_to_json_file();
    out << "\t\"iterations\": [\n";
    if (verbose) std::cout<<vtime()<<"Performing initial likelihood computation..."<<std::endl;
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
    //std::gamma_distribution<double> rgamma;

    //Writing initial state to output file
    state_to_json_file(0);

    for (int i=0;i<nbit;i++) {//Main MCMC loop
        if (!verbose && nbit > 0 && i > 0 && (100 * i) / nbit != (100 * (i - 1)) / nbit) print_progress(i, nbit);
        if (verbose) std::cout<<vtime()<<"Starting iteration "<<i<<std::endl;

        //#Adaptive MH move for rho using Gamma prior
        for (int rep=0;rep<updates[0];rep++) {
            rnorm.param(std::normal_distribution<double>::param_type(rho,sd[0]));
            double rho2=std::fabs(rnorm(rng));
            double prior_rho2=dgamma_mean_sd(rho2,prior_params[0],prior_params[1]);
            double ll2=tree->likelihood_anc_rec(rho2,delta,nu);
            double mh=ll2-ll+prior_rho2-prior_rho;
            if (verbose) std::cout<<vtime()<<std::scientific << std::setprecision(6)<<"Rho move proposed rho="<<rho<<" rho2="<<rho2<<" ll="<<ll<<" ll2="<<ll2<<" sd_rho="<<sd[0];
            if (log(runif(rng))<mh) {
                naccepted[0]++;
                if (verbose) std::cout<<" accepted\n";
                ll=ll2;
                rho=rho2;
                prior_rho=prior_rho2;
            } else {
                if (verbose) std::cout<<" rejected\n";
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
            if (verbose) std::cout<<vtime()<<std::scientific << std::setprecision(6)<<"Delta move proposed delta="<<delta<<" delta2="<<delta2<<" ll="<<ll<<" ll2="<<ll2<<" sd_delta="<<sd[1];
            if (log(runif(rng))<mh) {
                naccepted[1]++;
                if (verbose) std::cout<<" accepted\n";
                ll=ll2;
                delta=delta2;
                prior_delta=prior_delta2;
            } else {
                if (verbose) std::cout<<" rejected\n";
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
            if (verbose) std::cout<<vtime()<<std::scientific << std::setprecision(6)<<"Nu move proposed nu="<<nu<<" nu2="<<nu2<<" ll="<<ll<<" ll2="<<ll2<<" sd_nu="<<sd[2];
            if (log(runif(rng))<mh) {
                naccepted[2]++;
                if (verbose) std::cout<<" accepted\n";
                ll=ll2;
                nu=nu2;
                prior_nu=prior_nu2;
            } else {
                if (verbose) std::cout<<" rejected\n";
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
                lqr=lqr+std::log(pro/cur);
                /*double pro=2.0*std::sqrt(cur)+rnorm(rng);pro=pro*pro/4.0;
                lqr=lqr+0.5*std::log(pro/cur);*/
                /*double var=sd[3]*sd[3];
                rgamma.param(std::gamma_distribution<double>::param_type(cur/var,var));
                double pro=rgamma(rng);
                lqr=lqr+dgamma_shape_scale(cur,pro/var,var)-dgamma_shape_scale(pro,cur/var,var);*/
                tree2->get_node(b)->set_dist_father(pro);
            }
            tree2->jointml();
            double ll2=tree2->likelihood_anc_rec(rho,delta,nu);
            double prior_tree2=tree2->prior(this);
            double mh=ll2-ll+prior_tree2-prior_tree+lqr;
            if (verbose) std::cout<<vtime()<<std::scientific << std::setprecision(6)<<"Branch lengths move proposed ll="<<ll<<" ll2="<<ll2<<" lp="<<prior_tree<<" lp2="<<prior_tree2<<" mh="<<mh;
            if (log(runif(rng))<mh) {
                naccepted[3]++;
                if (verbose) std::cout<<" accepted\n";
                ll=ll2;
                prior_tree=prior_tree2;
                delete tree;
                tree=tree2;
            } else {
                if (verbose) std::cout<<" rejected\n";
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
            if (verbose) std::cout<<vtime()<<std::scientific << std::setprecision(6)<<"Tree move proposed ll="<<ll<<" ll2="<<ll2<<" lp="<<prior_tree<<" lp2="<<prior_tree2<<" mh="<<mh;
            if (log(runif(rng))<mh) {
                naccepted[4]++;
                if (verbose) std::cout<<" accepted\n";
                ll=ll2;
                prior_tree=prior_tree2;
                delete tree;
                tree=tree2;
            } else {
                if (verbose) std::cout<<" rejected\n";
                delete tree2;
            }
        }

        //Writing to output file
        if ((i+1)%thin==0) {
            out << "\t,\n";
            state_to_json_file(i+1);
        }

    }//End of main MCMC loop
    if (verbose==false && nbit>0) print_progress(nbit,nbit);
    out << "\t]\n";
    if (profiling) Profiler::prof.stop(__func__);
    if (verbose) std::cout<<vtime()<<"Finished MCMC."<<std::endl;
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
    out << "\t\t{\n"
        << std::scientific << std::setprecision(6)
        << "\t\t\t\"likelihood\": "        << ll                    << ",\n"
        << "\t\t\t\"prior_tree\": "        << prior_tree            << ",\n"
        << "\t\t\t\"prior_rho\": "         << prior_rho             << ",\n"
        << "\t\t\t\"prior_delta\": "       << prior_delta           << ",\n"
        << "\t\t\t\"prior_nu\": "          << prior_nu              << ",\n"
        << "\t\t\t\"sum_bra_len\": "       << tree->sum_bra_len()   << ",\n"
        << "\t\t\t\"tree\": \""           << tree->to_string()      << "\",\n"
        << "\t\t\t\"rho\": "               << rho                   << ",\n"
        << "\t\t\t\"delta\": "             << delta                 << ",\n"
        << "\t\t\t\"nu\": "                << nu                    << ",\n"
        << "\t\t\t\"sd_rho\": "            << sd[0]                 << ",\n"
        << "\t\t\t\"sd_delta\": "          << sd[1]                 << ",\n"
        << "\t\t\t\"sd_nu\": "             << sd[2]                 << ",\n"
        << "\t\t\t\"sd_lengths\": "        << sd[3]                 << ",\n"
        << "\t\t\t\"paccepted_rho\": "     << (iter * updates[0] == 0 ? 0.0 : 1.0 * naccepted[0] / (iter * updates[0])) << ",\n"
        << "\t\t\t\"paccepted_delta\": "   << (iter * updates[1] == 0 ? 0.0 : 1.0 * naccepted[1] / (iter * updates[1])) << ",\n"
        << "\t\t\t\"paccepted_nu\": "      << (iter * updates[2] == 0 ? 0.0 : 1.0 * naccepted[2] / (iter * updates[2])) << ",\n"
        << "\t\t\t\"paccepted_lengths\": " << (iter * updates[3] == 0 ? 0.0 : 1.0 * naccepted[3] / (iter * updates[3])) << ",\n"
        << "\t\t\t\"paccepted_spr\": "     << (iter * updates[4] == 0 ? 0.0 : 1.0 * naccepted[4] / (iter * updates[4])) << '\n'
        << "\t\t}\n";
}
