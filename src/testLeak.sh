#This is a bash script to test for memory leak using ASan, bypassing all CMake stuff
echo "constexpr const char* VERSION_STRING = \"0.0.1\";" > version.h
if command -v clang++ >/dev/null 2>&1; then
clang++ -g -O0 -fsanitize=address *.cpp -o /tmp/ClonalFrame2
else
g++ -g -O0 -fsanitize=address *.cpp -o /tmp/ClonalFrame2
fi
rm version.h
MallocNanoZone=0 ASAN_OPTIONS=detect_leaks=1 /tmp/ClonalFrame2 -v
MallocNanoZone=0 ASAN_OPTIONS=detect_leaks=1 /tmp/ClonalFrame2 -v ../data/simu.fasta /tmp/out
