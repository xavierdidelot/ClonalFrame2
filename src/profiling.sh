#This is a bash script to profile using gperftools and pprof, bypassing all CMake stuff
#Compilation happens only if there is no argument
if [ "$#" -eq 0 ]; then
  g++ -fopenmp -std=c++23 -O3 -g -fno-omit-frame-pointer -L/usr/local/lib -lprofiler *.cpp -o /tmp/ClonalFrame2profiling
fi 
CPUPROFILE=/tmp/profile.out /tmp/ClonalFrame2profiling -s 1 -S 10,1000,1,100,0.01,1e-4,1e-4 -y 2000
#pprof --pdf /tmp/ClonalFrame2profiling /tmp/profile.out > /tmp/profile.pdf
pprof -http=:8080 /tmp/ClonalFrame2profiling /tmp/profile.out

