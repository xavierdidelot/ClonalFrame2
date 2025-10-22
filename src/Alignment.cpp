//
// Created by Didelot, Xavier on 19/09/2025.
//

#include "Alignment.h"

#include <iostream>
#include <regex>

Alignment::~Alignment() {

}

Alignment::Alignment(const string &filename) {
    string line;
    ifstream file;
    if (verbose) cout<<"Opening file "<<filename<<endl;
    file.open(filename);
    if (file.is_open()==false) {
        cerr<<"Error opening file "<<filename<<endl;
        exit(1);
    }
    n=0;
    l=0;
    while (true) {
        getline(file,line);
        if (file.eof()) break;
        if (line.empty()) continue;
        if (line.at(0)=='#') continue;
        if (line.at(0)=='>') {n++;continue;}
        if (n==1) l=l+line.length();
    }
    file.clear();
    file.seekg(0);
    data=vector<vector<uint8_t>*>(n);
    for (size_t i=0;i<n;i++) data.at(i)=new vector<uint8_t>(l);
    names=vector<string>(n);
    n=0;
    l=0;
    regex reg("[> ]");
    while (true) {
        getline(file,line);
        if (file.eof()) break;
        if (line.empty()) continue;
        if (line.at(0)=='#') continue;
        if (line.at(0)=='>') {l=0;n++;names[n-1]=regex_replace(line,reg,"");continue;}
        for (size_t i=0;i<line.length();i++) {
            data.at(n-1)->at(l)=encode(line[i]);
            l++;
        }
    }
    if (verbose) cout<<"Finished reading file made of "<<n<<" genomes of length "<<l<<endl;
    file.close();
}

Alignment::Alignment(size_t num, size_t len) {
    n=num;
    l=len;
    data=vector<vector<uint8_t>*>(n);
    for (size_t i=0;i<n;i++) data.at(i)=new vector<uint8_t>(l,1);
}

Alignment::Alignment(Tree *t) {
    n=t->get_n();
    names=vector<string>(n);
    data=vector<vector<uint8_t>*>(n);
    for (size_t i=0;i<n;i++) {
        names[i]=t->get_node(i)->get_name();
        data[i]=t->get_node(i)->get_seq();
    }
    l=data[0]->size();
}

char Alignment::get_data(size_t i, size_t j) {return data[i]->at(j);}

size_t Alignment::get_n() {return n;}

size_t Alignment::get_l() {return l;}

uint8_t Alignment::encode(char in) {
    switch (in) {
        case 'a': return nuc_from_int(0);
        case 'A': return nuc_from_int(0);
        case 'c': return nuc_from_int(1);
        case 'C': return nuc_from_int(1);
        case 'g': return nuc_from_int(2);
        case 'G': return nuc_from_int(2);
        case 't': return nuc_from_int(3);
        case 'T': return nuc_from_int(3);
        default:  return nuc_from_int(4);
    }
}

char Alignment::decode(uint8_t in) {
    switch (in) {
        case 1u: return 'A';
        case 2u: return 'C';
        case 4u: return 'G';
        case 8u: return 'T';
        default: return 'N';
    }
}


void Alignment::write(const string &filename) {
    ofstream file;
    file.open(filename);
    if (file.is_open()==false) {
        cerr<<"Error opening file "<<filename<<endl;
        return;
    }
    for (size_t i=0;i<n;i++) {
        file<<">"<<names[i]<<endl;
        for (size_t j=0;j<l;j++) {
            file<<decode(get_data(i,j));
        }
        file<<endl;
    }
    file.close();
}

string Alignment::get_name(size_t i) {return names[i];}

double Alignment::dist(size_t i, size_t j) {
    int num=0,den=0;
    for (size_t a=0;a<l;a++) {
        if (nuc_is_missing(get_data(i,a))) continue;
        if (nuc_is_missing(get_data(j,a))) continue;
        den++;
        if (get_data(i,a)!=get_data(j,a)) num++;
    }
    return(-0.75*log(1.0-4.0/3.0*num/den));
}

vector<uint8_t> * Alignment::get_seq(size_t i) {
    return(data[i]);
}

void Alignment::miss_at_random(double x) {
    uniform_real_distribution unif(0.0, 1.0);
    for (size_t i=0;i<n;i++) {
        for (size_t j=0;j<l;j++) {
            if (unif(rng)<x) data[i]->at(j)=data[i]->at(j)|0b00010000;
        }
    }
}
