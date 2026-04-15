#This is a bash script to profile using gperftools and pprof, bypassing all CMake stuff
#Compilation happens only if there is no argument
if [ "$#" -eq 0 ]; then
  clang++ -std=c++23 -O2 -g -fno-omit-frame-pointer -L/usr/local/lib -lprofiler *.cpp -o /tmp/ClonalFrame2profiling
fi 
CPUPROFILE=/tmp/profile.out /tmp/ClonalFrame2profiling -s 1 -S -x 1000 -y 2000 -z 1
#pprof --pdf /tmp/ClonalFrame2 /tmp/profile.out > /tmp/profile.pdf
pprof -http=:8080 /tmp/ClonalFrame2profiling /tmp/profile.out

