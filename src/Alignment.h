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
    inline char getData(int i, int j) {return data[i][j];}
    inline int getN() {return N;}
    inline int getL() {return L;}
    void write(string filename);
    string getName(int i) {return names[i];}

private:
    int N;
    int L;
    vector<string> data;
    vector<string> names;
};


#endif //CF2_ALIGNMENT_H