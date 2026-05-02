#This is a bash script to do timing, bypassing all CMake stuff
#Compilation happens only if there is no argument
if [ "$#" -eq 0 ]; then
  g++-15 -fopenmp -std=c++23 -O3 *.cpp -o /tmp/ClonalFrame2
fi
/tmp/ClonalFrame2 -s 1 -S 10,1000,1e-4,1,100,0.01 -x 1000 -y 2000 -z 1

