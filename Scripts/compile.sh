#!/bin/bash

source ~/.bash_interface

#*# ************************************************************************ #*#
printHeader "Modifying sources..."

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

sed -e "s/%compiler%/${compiler}/" ./Sources/main.sed | sed -e "s/%target%/${target}/" | sed -e "s/%timestamp%/${timestamp}/" | sed -e "s/%platform%/${platform}/" > ./Sources/main.cc

rm -f ./Sources/main.sed

showProgressBar 100

#*# ************************************************************************ #*#
printHeader "Compiling..."

rm -f ./planar3D.out

#clang++ -Xpreprocessor -fopenmp -lomp -L/usr/local/opt/libomp/lib -I/usr/local/opt/libomp/include -std=c++11 -march=native -O3 -Weffc++ -Wall -pedantic -o ./planar3D.out ./Sources/*.cc

clang++ -std=c++11 -march=native -O3 -Weffc++ -Wall -pedantic -o ./planar3D.out ./Sources/*.cc

if [ ! -f ./planar3D.out ]; then
    printRed "Compilation failed."
    exit 1
fi

showProgressBar 100

#*# ************************************************************************ #*#
printHeader "Archiving..."

mkdir -p ./Archives

zip -rq ./Archives/$(date +'%d.%m.%Y_%H.%M.%S').zip ./ -i ./planar3D.out ./Sources/*.cc ./Sources/*.hh

rm -rf ./Sources

showProgressBar 100
