#This is a bash script to profile using gperftools and pprof, bypassing all CMake stuff
clang++ -O2 -g -fno-omit-frame-pointer -L/usr/local/lib -lprofiler *.cpp -o /tmp/ClonalFrame2
CPUPROFILE=/tmp/profile.out /tmp/ClonalFrame2 -v -S
pprof --pdf /tmp/ClonalFrame2 /tmp/profile.out > /tmp/profile.pdf
pprof -http=:8080 /tmp/ClonalFrame2 /tmp/profile.out

