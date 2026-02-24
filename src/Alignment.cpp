//
// Created by Didelot, Xavier on 19/09/2025.
//

#include "Alignment.h"

#include <iostream>
#include <regex>

using std::string,std::cerr,std::cout,std::endl,std::vector;

Alignment::~Alignment() {

}

Alignment * Alignment::read_fasta(const string &filename) {
    auto a = new Alignment();
    string line;
    std::ifstream file;
    if (verbose) cout<<"Opening file "<<filename<<endl;
    file.open(filename);
    if (file.is_open()==false) {
        cerr<<"Error opening file "<<filename<<endl;
        exit(1);
    }
    size_t n=0;
    size_t l=0;
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
    a->data=vector<vector<uint8_t>*>(n);
    for (size_t i=0;i<n;i++) a->data.at(i)=new vector<uint8_t>(l);
    a->names=vector<string>(n);
    n=0;
    l=0;
    std::regex reg("[> ]");
    while (true) {
        getline(file,line);
        if (file.eof()) break;
        if (line.empty()) continue;
        if (line.at(0)=='#') continue;
        if (line.at(0)=='>') {l=0;n++;a->names[n-1]=regex_replace(line,reg,"");continue;}
        for (size_t i=0;i<line.length();i++) {
            a->data.at(n-1)->at(l)=encode(line[i]);
            l++;
        }
    }
    if (verbose) cout<<"Finished reading file made of "<<n<<" genomes of length "<<l<<endl;
    file.close();
    a->n=n;
    a->l=l;
    a->pattern_analysis();
    return a;
}

Alignment::Alignment(size_t num, size_t len) {
    n=num;
    l=len;
    data=vector<vector<uint8_t>*>(n);
    for (size_t i=0;i<n;i++) data.at(i)=new vector<uint8_t>(l,1);
    pattern_analysis();
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
    pattern_analysis();
}

void Alignment::pattern_analysis() {
    if (verbose) cout<<"Starting pattern analysis"<<endl;
    auto data2=vector<vector<uint8_t>*>(n);
    for (size_t i=0;i<n;i++) data2.at(i)=new vector<uint8_t>(l);
    npat=0;
    patindex=vector<size_t>(l,0);
    weights =vector<size_t>(l,0);
    for (size_t j=0;j<l;j++) {
        bool toadd=true;
        size_t j2=0;
        for (j2=0;j2<npat;j2++) {
            bool samepat=true;
            for (size_t i=0;i<n;i++) {
                if (data[i]->at(j)!=data2[i]->at(j2)) {samepat=false;break;}
            }
            if (samepat) {toadd=false;break;}
        }
        if (toadd) {
            for (size_t i=0;i<n;i++) {
                data2[i]->at(npat)=data[i]->at(j);
            }
            j2=npat;
            npat++;
        }
        patindex[j]=j2;
        weights[j2]++;
    }
    weights.resize(npat);
    for (size_t i=0;i<n;i++) data2[i]->resize(npat);
    data=data2;
    if (verbose) cout<<"Finished pattern analysis and found "<<npat<<" patterns"<<endl;
}

char Alignment::get_data(size_t i, size_t j) {return data[i]->at(patindex[j]);}

size_t Alignment::get_n() {return n;}

size_t Alignment::get_l() {return l;}

size_t Alignment::get_npat() {return npat;}

uint8_t Alignment::encode(char in) {
    switch (in) {
        case 'a': return 0b00000001;
        case 'A': return 0b00000001;
        case 'c': return 0b00000010;
        case 'C': return 0b00000010;
        case 'g': return 0b00000100;
        case 'G': return 0b00000100;
        case 't': return 0b00001000;
        case 'T': return 0b00001000;
        default:  return 0b00010000;
    }
}

char Alignment::decode(uint8_t in) {
    switch (in) {
        case 0b00000001: return 'A';
        case 0b00000010: return 'C';
        case 0b00000100: return 'G';
        case 0b00001000: return 'T';
        default        : return 'N';
    }
}


void Alignment::write(const string &filename) {
    std::ofstream file;
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
    size_t num=0,den=0;
    for (size_t a=0;a<l;a++) {
        if (nuc_is_unsure(get_data(i,a))) continue;
        if (nuc_is_unsure(get_data(j,a))) continue;
        den+=weights[j];
        if (get_data(i,a)!=get_data(j,a)) num+=weights[j];
    }
    return(-0.75*log(1.0-4.0/3.0*num/den));
}

vector<uint8_t> * Alignment::get_seq(size_t i) {
    return(data[i]);
}

void Alignment::miss_at_random(double x) {
    std::uniform_real_distribution unif(0.0, 1.0);
    for (size_t i=0;i<n;i++) {
        for (size_t j=0;j<npat;j++) {
            if (unif(rng)<x) data[i]->at(j)=data[i]->at(j)|0b00010000;
        }
    }
}
