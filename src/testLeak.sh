#This is a bash script to test for memory leak using ASan, bypassing all CMake stuff
echo "constexpr const char* VERSION_STRING = \"0.0.1\";" > version.h
g++ -o /tmp/ClonalFrame2 -O0 -fsanitize=leak *.cpp
rm version.h
/tmp/ClonalFrame2
