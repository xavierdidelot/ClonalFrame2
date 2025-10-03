#This is a bash script to test for memory leak using ASan, bypassing all CMake stuff
echo "constexpr const char* VERSION_STRING = \"0.0.1\";" > version.h
g++ -g -O0 -fsanitize=address *.cpp -o /tmp/ClonalFrame2
rm version.h
/tmp/ClonalFrame2
/tmp/ClonalFrame2 ../data/simu.fasta /tmp/out
