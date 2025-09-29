#include <iostream>
#include <unistd.h>

#include "Move.h"
#include "Param.h"
#include "Alignment.h"

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
    srand(time(nullptr));
    auto * p=new Param();//TEST LEAK
    p->set_alpha(0);
    int * x=(int *)malloc(1000);
    x[1]=1;
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
        delete(a);
        return 0;
}