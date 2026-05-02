#This is a bash script to test for memory leak using ASan, bypassing all CMake stuff
#Compilation happens only if there is no argument
if [ "$#" -eq 0 ]; then
  clang++ -fopenmp -std=c++23 -g -O1 -fsanitize=address -fno-omit-frame-pointer *.cpp -o /tmp/ClonalFrame2leaking
fi
MallocNanoZone=0 ASAN_OPTIONS=detect_leaks=1 /tmp/ClonalFrame2leaking -S 10,1000,1e-4,1,100,0.01 -v -x 10 -y 20 -z 1
