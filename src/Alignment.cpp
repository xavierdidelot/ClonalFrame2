//
// Created by Didelot, Xavier on 19/09/2025.
//

#include "Alignment.h"

#include <iostream>

Alignment::~Alignment() {
    delete(data);
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
    int n=0;
    int l=0;
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
    cout<<n<<" "<<l<<endl;
    data=gsl_matrix_char_calloc(n,l);
    n=-1;
    l=0;
    while (true) {
        getline(file,line);
        if (file.eof()) break;
        if (line.empty()) continue;
        if (line.at(0)=='#') continue;
        if (line.at(0)=='>') {l=0;n++;continue;}
        for (int i=0;i<line.length();i++) {
            gsl_matrix_char_set(data,n,l++,standardise(line[i++]));
        }
    }
    file.close();
}

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

char Alignment::getData(int i,int j) {
    return(gsl_matrix_char_get(data,i,j));
}

int Alignment::getN() {
    return(data->size1);
}

int Alignment::getL() {
    return(data->size2);
}