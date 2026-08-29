# ClonalFrame2

## Introduction

ClonalFrame2 is a software to perform Bayesian phylogenetics 
accounting for recombination in whole bacterial genomes.
It is a new implementation of the original ClonalFrame
algorithm described in 
[Didelot and Falush (2007)](https://doi.org/10.1534/genetics.106.063305).

## Installation

First you need a local copy of this repository.
You can clone it using the command:
```
git clone https://github.com/xavierdidelot/ClonalFrame2
```
or you can download and unpack a zip file containing the 
[latest release](https://github.com/xavierdidelot/ClonalFrame2/archive/refs/heads/main.zip). 

The recommended method of compilation is to run the following
command in the directory where the repository was copied:
```
make
```

The command above uses `make` and `cmake` for configuration. 
Alternatively compilation can be achieved directly with:
```
g++ -fopenmp -std=c++17 -O3 src/*.cpp -o ClonalFrame2
```

If OpenMP is not available, ClonalFrame2 can be compiled
for non-parallel usage with:
```
g++ -std=c++17 -O3 src/*.cpp -o ClonalFrame2
```

Compilation will generate an executable called
`ClonalFrame2`. You can make it available system-wide
using:
```
make install
```
or by copying it to a location included in your `$PATH`.

## Usage from the command line

The `ClonalFrame2` executable has the following options:

```
    -i FILE     Select the input file
    -o FILE     Select the output file
    -x NUM      Sets the total number of iterations (default is 10000)
    -y NUM      Sets the number of iterations between samples (default is x/2000)
    -s NUM      Use given seed to initiate random number generator
    -c NUM      Number of cores to use (default is all available)
    -S STR      Generate a simulated dataset and run on it
                STR is 7 comma-separated values for number of sequences,
                length of sequences, rho, delta, nu, branch length mean and sd
    -I STR      Specify starting parameter values (default is 1,500,0.05)
                STR is 3 comma-separated values for rho, delta and nu
                If -S is used without -I then start from correct values and tree
    -U STR      Indicates number of updates to use per iteration for each move
                STR is 5 int values for rho, delta, nu, lengths and SPR
    -P STR      Indicates prior as 8 comma-separated values for mean and sd
                of rho, delta, nu and b (default 1,1,500,500,0.05,0.05,1e-4,1e-4)
    -p          Activate profiler
    -v          Verbose mode
    -h          Show help
```

For example, to run on a simulated dataset with 100 sequences of length 1Mbp, with
parameters rho=0.1, delta=1000bp, nu=0.01 and branch lengths
exponentially distributed with mean 1e-4:

```
ClonalFrame2 -S 100,1000000,0.1,1000,0.01,1e-4,1e-4 -o out.json
```

For example, to run on a FASTA file called `data.fasta`:

```
ClonalFrame2 -i data.fasta -o out.json
```

All results are stored in the json output file. The format
is relatively straightforward, but most users
will want to analyse results using the R package `RClonalFrame2` 
described below.

## The RClonalFrame2 package

ClonalFrame2 comes with an R package called `RClonalFrame2` which
provides functions to run from R and analyse results.
You need to install `ClonalFrame2` as described above
before you install the R package.
You can then install `RClonalFrame2` using the following command from the directory
where you copied the repository:

```
R CMD INSTALL RClonalFrame2
```

Once installed, you can load it in R using:

```
library(RClonalFrame2)
```

You can then run ClonalFrame2 from R, for example to run
on the same simulation conditions as previously described you
can use: 
```
x=run_cf2(simu_pars=c(100,1000000,0.1,1000,0.01,1e-4,1e-4),outputfile='out.json')
```
or to run on a FASTA file called `data.fasta` you can use:
```
x=run_cf2(inputfile='data.fasta',outputfile='out.json')
```

You can also load the results of a previous run using:
```
x=load_cf2('out.json')
```

The results can then be analysed and plotted using several functions,
for example:
```
print(x)                #shows parameters
plot(x)                 #plots parameters
plot(x,'trace')         #plots traces
m=medoid(x)             #return index of medoid tree
x$iterations$tree[m]    #newick string of medoid tree
plot(x,'recombination') #plots recombination events on medoid tree
```

For more details about the functionality of the RClonalFrame2 package
see the help using:
```
help("RClonalFrame2")
```

## Getting help

If you have any problem or question please create an issue
[here](https://github.com/xavierdidelot/ClonalFrame2/issues) or get in
touch by emailing `xavier.didelot@gmail.com`