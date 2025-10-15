//
// Created by Didelot, Xavier on 19/09/2025.
//

#include "Alignment.h"

#include <iostream>
#include <regex>

Alignment::~Alignment() {
    //for (int i = 0; i < n; i++) {delete data->at(i);}//An alignment does not delete the seqs: this is done by the tree
    delete data;
}

Alignment::Alignment(string filename) {
    string line;
    ifstream file;
    cout<<"Opening file "<<filename<<endl;
    file.open(filename);
    if (file.is_open()==false) {
        cout<<"Error opening file "<<filename<<endl;
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
    data=new vector<vector<uint8_t>*>(n);
    for (int i=0;i<n;i++) data->at(i)=new vector<uint8_t>(l);
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
        for (unsigned int i=0;i<line.length();i++) {
            data->at(n-1)->at(l)=encode(line[i]);
            l++;
        }
    }
    cout<<"Finished reading file made of "<<n<<" genomes of length "<<l<<endl;
    file.close();
}

Alignment::Alignment(int num, int len) {
    n=num;
    l=len;
    data=new vector<vector<uint8_t>*>(n);
    for (int i=0;i<n;i++) data->at(i)=new vector<uint8_t>(l,1);
}

Alignment::Alignment(Tree *t) {
    n=t->get_n();
    names=vector<string>(n);
    data=new vector<vector<uint8_t>*>(n);
    for (int i=0;i<n;i++) {
        names[i]=t->get_node(i)->get_name();
        data->at(i)=t->get_node(i)->get_seq();
    }
    l=data->at(0)->size();
}

char Alignment::get_data(int i, int j) {return data->at(i)->at(j);}

int Alignment::get_n() {return n;}

int Alignment::get_l() {return l;}

uint8_t Alignment::encode(char in) {
    switch (in) {
        case 'a':
            return 1u;
        case 'A':
            return 1u;
        case 'c':
            return 1u << 1;
        case 'C':
            return 1u << 1;
        case 'g':
            return 1u << 2;
        case 'G':
            return 1u << 2;
        case 't':
            return 1u << 3;
        case 'T':
            return 1u << 3;
        default:
            return 1u << 4;
    }
}

char Alignment::decode(uint8_t in) {
    switch (in) {
        case 1u:
            return 'A';
        case 1u << 1:
            return 'C';
        case 1u << 2:
            return 'G';
        case 1u << 3:
            return 'T';
        default:
            return 'N';
    }
}


void Alignment::write(string filename) {
    ofstream file;
    file.open(filename);
    if (file.is_open()==false) {
        cout<<"Error opening file "<<filename<<endl;
        exit(1);
    }
    for (int i=0;i<n;i++) {
        file<<">"<<names[i]<<endl;
        for (int j=0;j<l;j++) {
            file<<decode(get_data(i,j));
        }
        file<<endl;
    }
    file.close();
}

string Alignment::get_name(int i) {return names[i];}

double Alignment::dist(int i, int j) {
    int num=0,den=0;
    for (int a=0;a<l;a++) {
        if (nuc_is_missing(get_data(i,a))) continue;
        if (nuc_is_missing(get_data(j,a))) continue;
        den++;
        if (get_data(i,a)!=get_data(j,a)) num++;
    }
    return(-0.75*log(1.0-4.0/3.0*num/den));
}

vector<uint8_t> * Alignment::get_seq(int i) {
    return(data->at(i));
}
