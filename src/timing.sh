#This is a bash script to do timing, bypassing all CMake stuff
#If there is no argument then compilation happens and 12 cores are used
#If there is an argument then it is the number of cores and no compilation happens
if [ "$#" -eq 0 ]; then
  g++ -fopenmp -std=c++23 -O3 *.cpp -o /tmp/ClonalFrame2
  ncores=12
else
  ncores=$1
fi
/tmp/ClonalFrame2 -s 1 -S 100,1000000,1e-4,1,100,0.01 -x 1000 -y 2000 -z 1 -c $ncores -p

