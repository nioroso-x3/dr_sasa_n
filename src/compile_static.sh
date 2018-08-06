#!/bin/sh
#linux
g++ -static -static-libgcc -static-libstdc++ -O3 -o ../build/dr_sasa.bin -std=c++11 *.cpp
#windows64
x86_64-w64-mingw32-g++ -static -static-libgcc -static-libstdc++ -O3 -o ../build/dr_sasa.exe -std=c++11 *.cpp
