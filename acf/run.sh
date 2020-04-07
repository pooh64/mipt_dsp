g++ -std=c++14 -O2 -o task task.cpp
mkdir -p out
mkdir -p png
./task
gnuplot plotscript
