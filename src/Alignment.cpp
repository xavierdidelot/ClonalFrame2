//
// Created by Didelot, Xavier on 19/09/2025.
//

#include "Alignment.h"

#include <iostream>
#include <regex>

Alignment::~Alignment() {

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
    this->n=n;
    this->l=l;
    data=vector<string>(n);
    for (int i=0;i<n;i++) {string s(l,'*');data.push_back(s);}
    names=vector<string>(n);
    n=-1;
    l=0;
    regex reg("[> ]");
    while (true) {
        getline(file,line);
        if (file.eof()) break;
        if (line.empty()) continue;
        if (line.at(0)=='#') continue;
        if (line.at(0)=='>') {l=0;n++;names[n]=regex_replace(line,reg,"");continue;}
        for (int i=0;i<line.length();i++) {
            data[n][l++]=standardise(line[i]);
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
}
