#This is a bash script to do timing, bypassing all CMake stuff
#Compilation happens only if there is no argument
if [ "$#" -eq 0 ]; then
  clang++ -std=c++23 -O3 *.cpp -o /tmp/ClonalFrame2
fi
/tmp/ClonalFrame2 -s 1 -S -x 1000 -y 2000 -z 1

