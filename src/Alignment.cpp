//
// Created by Didelot, Xavier on 19/09/2025.
//

#include "Alignment.h"

#include <iostream>
#include <regex>

Alignment::~Alignment() {
    for (size_t i=0;i<data->size();i++) delete data->at(i);
    delete data;
    delete names;
    delete patindex;
    delete weights;
    delete refsites;
    delete refsites_diff;
}

Alignment * Alignment::read_fasta(const string &filename) {
    auto a = new Alignment();
    string line;
    std::ifstream file;
    if (verbose) vprint("Opening file {}",filename);
    file.open(filename);
    if (file.is_open()==false) {
        std::println("Error opening file {}",filename);
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
    a->data=new vector<vector<uint8_t>*>(n);
    for (size_t i=0;i<n;i++) a->data->at(i)=new vector<uint8_t>(l);
    a->names=new vector<string>(n);
    n=0;
    l=0;
    std::regex reg("[> ]");
    while (true) {
        getline(file,line);
        if (file.eof()) break;
        if (line.empty()) continue;
        if (line.at(0)=='#') continue;
        if (line.at(0)=='>') {l=0;n++;a->names->at(n-1)=regex_replace(line,reg,"");continue;}
        for (size_t i=0;i<line.length();i++) {
            a->data->at(n-1)->at(l)=encode(line[i]);
            l++;
        }
    }
    if (verbose) vprint("Finished reading file made of {} genomes of length {}",n,l);
    file.close();
    a->n=n;
    a->l=l;
    a->pattern_analysis();
    return a;
}

Alignment* Alignment::make_nonpoly(size_t n, size_t l) {
    auto a = new Alignment();
    a->n=n;
    a->l=l;
    a->data=new vector<vector<uint8_t>*>(n);
    for (size_t i=0;i<n;i++) a->data->at(i)=new vector<uint8_t>(l,1);
    a->pattern_analysis();
    return a;
}

Alignment* Alignment::make_from_tree(Tree *t, bool include_ancestors) {
    auto a = new Alignment();
    size_t n=t->get_n();
    if (include_ancestors) {
        n=n*2-1;
        for (size_t i=t->get_n();i<n;i++) t->get_node(i)->set_name("NODE"+std::to_string(i));
    }
    a->n=n;
    a->names=new vector<string>(n);
    a->data=new vector<vector<uint8_t>*>(n);
    for (size_t i=0;i<n;i++) {
        a->names->at(i)=t->get_node(i)->get_name();
        a->data->at(i)=new vector(*t->get_node(i)->get_seq());
    }
    size_t l=a->data->at(0)->size();
    if (t->get_alignment()==nullptr) {
        a->l=l;
        a->npat=l;
        a->weights=new vector<size_t>(l,1);
        a->patindex=new vector<size_t>(l,0);
        for (size_t j=0;j<l;j++) a->patindex->at(j)=j;
    } else {
        a->l=t->get_alignment()->get_l();
        a->npat=t->get_alignment()->get_npat();
        a->weights=new vector(*t->get_alignment()->get_weights());
        a->patindex=new vector(*t->get_alignment()->get_patindex());
    }
    a->make_refsites();
    return a;
}

void Alignment::pattern_analysis() {
    if (verbose) vprint("Starting pattern analysis");
    auto data2=new vector<vector<uint8_t>*>(n);
    size_t lplus=l;
    if (lplus<100) lplus=100;
    for (size_t i=0;i<n;i++) data2->at(i)=new vector<uint8_t>(lplus);
    npat=0;
    delete patindex;
    delete weights;
    patindex=new vector<size_t>(l,0);
    weights =new vector<size_t>(lplus,0);
    std::unordered_map<size_t, vector<size_t>> hashmap;
    for (int j=0;j<4;j++) {//Adding the four non-polymorphic patterns as the first four patterns
        size_t h = 0;
        for (size_t i = 0; i < n; i++)
            h = h * 131 + (1u<<j);
        for (size_t i=0;i<n;i++)
            data2->at(i)->at(npat)=1u<<j;
        hashmap[h].push_back(npat);
        npat++;
    }

    for (size_t j=0;j<l;j++) {
        // Compute hash
        size_t h = 0;
        for (size_t i = 0; i < n; i++)
            h = h * 131 + data->at(i)->at(j);
        bool toadd=true;
        size_t j2=0;
        // Check candidates with same hash (collision-safe)
        auto it = hashmap.find(h);
        if (it != hashmap.end()) {
            for (size_t cand : it->second) {
                bool same = true;
                for (size_t i = 0; i < n; i++)
                    if (data->at(i)->at(j) != data2->at(i)->at(cand)) { same = false; break; }
                if (same) { toadd = false; j2 = cand; break; }
            }
        }        if (toadd) {
            for (size_t i=0;i<n;i++) {
                data2->at(i)->at(npat)=data->at(i)->at(j);
            }
            hashmap[h].push_back(npat);
            j2=npat;
            npat++;
        }
        patindex->at(j)=j2;
        weights->at(j2)++;
    }
    weights->resize(npat);
    for (size_t i=0;i<n;i++) data2->at(i)->resize(npat);
    for (size_t i=0;i<n;i++) delete data->at(i);
    delete data;
    data=data2;
    make_refsites();
    if (verbose) vprint("Finished pattern analysis and found {} patterns",npat);
}

void Alignment::make_refsites() {
    delete refsites;
    delete refsites_diff;
    refsites=new vector<size_t>(0,0);
    refsites_diff=new vector<size_t>(0,0);
    int dist=0;
    for (size_t i=0;i<l;i++) {
        if (dist==maxskip-1 || patindex->at(i)>3 || i==l-1) {
            refsites->push_back(i);
            refsites_diff->push_back(dist);
            dist=0;
        } else {
            dist++;
        }
    }
}

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


void Alignment::write(const string &filename, bool include_uncertainty) {
    std::ofstream file;
    file.open(filename);
    if (file.is_open()==false) {
        std::println("Error opening file {}",filename);
        return;
    }
    for (size_t i=0;i<n;i++) {
        std::println(file,">{}",names->at(i));
        for (size_t j=0;j<l;j++) {
            if (include_uncertainty) std::print(file,"{}",decode(get_data(i,j)));
            else std::print(file,"{}",decode(get_data(i,j)&0b00001111));
        }
        std::println(file,"");
    }
    file.close();
}

double Alignment::dist(size_t i, size_t j) {
    size_t num=0,den=0;
    for (size_t a=0;a<npat;a++) {
        if (nuc_is_unsure(data->at(i)->at(a))) continue;
        if (nuc_is_unsure(data->at(j)->at(a))) continue;
        den+=weights->at(a);
        if (data->at(i)->at(a)!=data->at(j)->at(a)) num+=weights->at(a);
    }
    return(-0.75*log(1.0-4.0/3.0*num/den));
}

void Alignment::miss_at_random(double x) {
    std::uniform_real_distribution unif(0.0, 1.0);
    for (size_t i=0;i<n;i++) {
        for (size_t j=0;j<npat;j++) {
            if (unif(rng)<x) data->at(i)->at(j)=data->at(i)->at(j)|0b00010000;
        }
    }
}

vector<vector<double>> * Alignment::distmat() {
    auto d=new vector<vector<double>>(n, vector<double>(n, 0.0));
    for (size_t i = 0; i < n; i++) for (size_t j = 0; j < i; j++) {
        (*d)[i][j] = dist(i, j);
        (*d)[j][i] = (*d)[i][j];
    }
    return d;
}
