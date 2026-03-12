//
// Created by Didelot, Xavier on 21/10/2025.
//

#include "Param.h"

using std::string,std::cerr,std::cout,std::endl,std::vector;

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

    cout << "\r[" << bar << "] " << static_cast<int>(fraction * 100) << "%  " << std::flush;

    if (current == total) cout << '\n';
}

void Param::mcmc(int burn, int nbit, int thin, string outfile) {
    std::ofstream file;
    if (!outfile.empty()) {
        file.open(outfile);
        if (file.is_open() == false) cerr << "Error opening file " << outfile << endl;
    }

    for (int i=0;i<nbit;i++) {
        if (verbose==false && floor(i/nbit*100)==0) print_progress(i,nbit);
        if (verbose) cout<<"Starting iteration "<<i<<endl;
        ll=tree->likelihood();
        tree->jointml();
        ll=tree->likelihood_anc();
        if (i>burn && i%thin==0) {
            file << tree->to_string() << endl;
        }
    }
    if (verbose==false) print_progress(nbit,nbit);
    cout<<endl;
    file.close();
}
