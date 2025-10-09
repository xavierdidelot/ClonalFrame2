//
// Created by Didelot, Xavier on 19/09/2025.
//

#include "Alignment.h"

#include <iostream>
#include <regex>

const char * lets="ACGT";

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
    data=new vector<string*>(n);
    for (int i=0;i<n;i++) data->at(i)=new string(l,'*');
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
            data->at(n-1)->at(l)=standardise(line[i]);
            l++;
        }
    }
    cout<<"Finished reading file made of "<<n<<" genomes of length "<<l<<endl;
    file.close();
}

Alignment::Alignment(int num, int len) {
    n=num;
    l=len;
    data=new vector<string*>(n);
    for (int i=0;i<n;i++) data->at(i)=new string(l,'A');
}

Alignment::Alignment(Tree *t) {
    n=t->get_n();
    names=vector<string>(n);
    data=new vector<string*>(n);
    for (int i=0;i<n;i++) {
        names[i]=t->get_node(i)->get_name();
        data->at(i)=t->get_node(i)->get_seq();
    }
    l=data->at(0)->length();
}

char Alignment::get_data(int i, int j) {return data->at(i)->at(j);}

int Alignment::get_n() {return n;}

int Alignment::get_l() {return l;}

char Alignment::standardise(char in) {
    switch (in) {
        case 'a':
            return 'A';
        case 'A':
            return 'A';
        case 't':
            return 'T';
        case 'T':
            return 'T';
        case 'c':
            return 'C';
        case 'C':
            return 'C';
        case 'g':
            return 'G';
        case 'G':
            return 'G';
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
            file<<get_data(i,j);
        }
        file<<endl;
    }
    file.close();
}

string Alignment::get_name(int i) {return names[i];}

double Alignment::dist(int i, int j) {
    double d=0;
    for (int a=0;a<l;a++)
        if (data->at(i)->at(a)!=data->at(j)->at(a)) d++;
    return(d/l);
}

string * Alignment::get_seq(int i) {
    return(data->at(i));
}
