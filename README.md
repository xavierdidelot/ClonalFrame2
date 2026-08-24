# ClonalFrame2

## Installation

In the directory where the source was downloaded, simply run:
```
make
```

## Usage

```
Usage: ClonalFrame2 [OPTIONS]

    Options:
    -i FILE     Select the input file
    -o FILE     Select the output file
    -x NUM      Sets the number of burn-in iterations (default is y/2)
    -y NUM      Sets the total number of iterations (default is 10000)
    -z NUM      Sets the number of iterations between samples (default is y/2000)
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
