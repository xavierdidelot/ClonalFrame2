#include "../src/Tree.h"
#include "../src/Alignment.h"
#include <cassert>
#include <iostream>

#include <string>

using namespace std;

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    auto t = Tree::random_tree(1000,0.0001);
    int m=t->simulate_sequences(10000);
    int s1=0,s2=0;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << "Simulation time: " << elapsed.count() << " seconds" << endl;
    cout<<"Simulated number of substitutions: "<<m<<endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i=0;i<100;i++) {
        s1=t->parsimony();
    }
    elapsed = std::chrono::high_resolution_clock::now() - start;
    cout << "Parsimony time: " << elapsed.count() << " seconds" << endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i=0;i<100;i++) {
        s2=t->parsimony_recursive();
    }
    elapsed = std::chrono::high_resolution_clock::now() - start;
    cout << "Parsimony_recursive time: " << elapsed.count() << " seconds" << endl;

    cout<<s1<<" "<<s2<<endl;
    delete(t);
}