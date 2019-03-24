#!/bin/bash

source ~/.bash_interface

printHeader "Compiling..."

timestamp=$(date +'%d.%m.%Y %H:%M:%S')

compiler="unknown"
target="unknown"
platform="unknown"

mkdir -p ./Sources

rm -f ./Sources/*.cc ./Sources/*.hh

cp ../Sources/Sequential/*.cc ./Sources/
cp ../Sources/Sequential/*.hh ./Sources/

cp /usr/local/include/ai.hh ./Sources/

mv ./Sources/main.cc ./Sources/main.sed

sed -e "s/%compiler%/${compiler}/g" ./Sources/main.sed | sed -e "s/%target%/${target}+debug/g" | sed -e "s/%timestamp%/${timestamp}/g" | sed -e "s/%platform%/${platform}/g" > ./Sources/main.cc

rm ./Sources/main.sed

clang++ -std=c++11 -O0 -g -Weffc++ -Wextra -Werror -pedantic -o ./planar3D.out ./Sources/*.cc
