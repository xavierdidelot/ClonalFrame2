#This is a bash script to test for memory leak using ASan, bypassing all CMake stuff
if command -v clang++ >/dev/null 2>&1; then
CXX=clang++
else
CXX=g++
fi
$CXX -g -O0 -fsanitize=address *.cpp -o /tmp/ClonalFrame2
MallocNanoZone=0 ASAN_OPTIONS=detect_leaks=1 /tmp/ClonalFrame2 -S -v
