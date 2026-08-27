# ClonalFrame2

## Compilation

The recommended method of compilation is to run the following
command in the directory where the repository was downloaded:
```
make
```

The command above uses make and cmake for configuration. 
Alternatively compilation can be achieved directly with:
```
g++ -fopenmp -std=c++17 -O3 src/*.cpp -o ClonalFrame2
```

If OpenMP is not available, ClonalFrame2 can be compiled
for non-parallel usage with:
```
g++ -std=c++17 -O3 src/*.cpp -o ClonalFrame2
```

## Usage

```
Usage: ClonalFrame2 [OPTIONS]

    Options:
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
