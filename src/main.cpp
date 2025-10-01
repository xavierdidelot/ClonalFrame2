#include <iostream>
#include <unistd.h>
#include <thread>
#include <chrono>

#include "Alignment.h"
#include "Tree.h"

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
    std::this_thread::sleep_for(std::chrono::seconds(1)); // Pause for profiling
    srand(time(nullptr));
    auto * tr=new Tree(10);//CREATE A LEAK
    cout<<tr->to_string()<<endl;
    int burn=10000,nbit=10000,thin=10;
    int c;
    while ((c = getopt (argc, argv, "x:y:z:s:v")) != -1)
        switch (c)
        {
        case('x'):burn=atoi(optarg);break;
        case('y'):nbit=atoi(optarg);break;
        case('z'):thin=atoi(optarg);break;
        case('s'):srand(atoi(optarg));break;
        case('v'):verbose=true;break;
        case '?':cout<<"Unknown option."<<endl;return 1;
        default:
                abort ();
        }

    if (argc-optind==0) {cout<<help<<endl;exit(0);}
    if (argc-optind<2) {cout<<"Too few non option arguments."<<endl;return 1;}
        if (argc-optind>2) {cout<<"Too many non option arguments."<<endl;return 1;}

        string inputfile=string(argv[optind++]);
        string outfile =string(argv[optind++]);
        auto * a = new Alignment(inputfile);
        auto * t = new Tree(a->get_n());
        cout << t->to_string() << endl;
        delete(a);
        return 0;
}