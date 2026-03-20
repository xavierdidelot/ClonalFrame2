#This is a bash script to test for memory leak using ASan, bypassing all CMake stuff
clang++ -std=c++23 -g -O1 -fsanitize=address -fno-omit-frame-pointer *.cpp -o /tmp/ClonalFrame2
MallocNanoZone=0 ASAN_OPTIONS=detect_leaks=1 /tmp/ClonalFrame2 -S -v -x 100 -y 200 -z 1
