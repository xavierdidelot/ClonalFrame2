#include <iostream>
#include <unistd.h>

#include "Alignment.h"
#include "Tree.h"
#include "RNGGlobal.h"
#include "version.h"

static string help=
    R"(Usage: ClonalFrame2 [OPTIONS] inputfile outputfile

    Options:
    -x NUM      Sets the number of burn-in iterations (default is 10000)
    -y NUM      Sets the number of iterations after burn-in (default is 10000)
    -z NUM      Sets the number of iterations between samples (default is 10)
    -s NUM      Use given seed to initial random number generator
    -v          Verbose mode
)";


bool verbose=false;

int main(int argc, char *argv[]) {
    cout<<"This is ClonalFrame2, version "<<VERSION_STRING<<endl;

    int* array = new int[100];array[0] = 42;//MEMORY LEAK
    int burn=10000,nbit=10000,thin=10;
    int c;
    while ((c = getopt (argc, argv, "x:y:z:s:v")) != -1)
        switch (c)
        {
        case('x'):burn=stoi(optarg);break;
        case('y'):nbit=stoi(optarg);break;
        case('z'):thin=stoi(optarg);break;
        case('s'):setSeed(stoi(optarg));break;
        case('v'):verbose=true;break;
        case '?':cout<<"Unknown option."<<endl;return 0;
        default:
                abort ();
        }

    if (argc-optind==0) {cout<<help<<endl;return 0;}
    if (argc-optind<2) {cout<<"Too few non option arguments."<<endl;return 0;}
        if (argc-optind>2) {cout<<"Too many non option arguments."<<endl;return 0;}

        string inputfile=string(argv[optind++]);
        string outfile =string(argv[optind++]);
        auto * a = new Alignment(inputfile);
        auto * t = new Tree(a->get_n());
        t -> nj(a);
        cout << t->to_string() << endl;
        delete(t);
        delete(a);
        cout << burn<< " "<<nbit<<" "<<thin<<" "<<outfile<<endl;
        return 0;
}