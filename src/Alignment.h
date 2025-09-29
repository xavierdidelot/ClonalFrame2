//
// Created by Didelot, Xavier on 19/09/2025.
//

#ifndef CF2_ALIGNMENT_H
#define CF2_ALIGNMENT_H

#include <string>
#include <fstream>
#include <vector>
using namespace std;

class Alignment {
public:
    ~Alignment();
    Alignment(string filename);
    char standardise(char in);
    inline char get_data(int i, int j) {return data[i][j];}
    inline int get_n() {return n;}
    inline int get_l() {return l;}
    void write(string filename);
    string get_name(int i) {return names[i];}

private:
    int n;
    int l;
    vector<string> data;
    vector<string> names;
};


#endif //CF2_ALIGNMENT_H