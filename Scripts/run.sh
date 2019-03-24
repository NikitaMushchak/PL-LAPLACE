#!/bin/bash

source ~/.bash_interface

#*# ************************************************************************ #*#
printHeader "Setting enviroment..."

mkdir -p ./InitialConditions

cp -f ../InitialConditions/*.txt ./InitialConditions/

mkdir -p ./Results

rm -rf ./Results/*

mkdir ./Results/Opening
mkdir ./Results/Pressure
mkdir ./Results/Flows

showProgressBar 100

#*# ************************************************************************ #*#
printHeader "Running..."

./planar3D.out --Q=0.08 --mu=0.4 --n=1 --E=32 --C=0 --time=10 --ts=60 --mesh-size=4 --cell-size=5 --barriers ${*}
