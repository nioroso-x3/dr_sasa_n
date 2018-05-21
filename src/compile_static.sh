#!/bin/sh
#linux
g++ -static -static-libgcc -static-libstdc++ -O3 -o dr_sasa -std=c++11 *.cpp
#windows
x86_64-w64-mingw32-g++ -static -static-libgcc -static-libstdc++ -O3 -o dr_sasa.exe -std=c++11 *.cpp
