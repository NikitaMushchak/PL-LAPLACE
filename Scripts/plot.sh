#!/bin/bash

files=(./Results/Opening/time*_m.txt)
start=1
j=3
filesCounter=$(ls -1 ./Results/Opening/time*_m.txt | wc -l)

gnuplot << EOD
reset
set terminal gif animate size 1200,1200 enhanced font 'Verdana,10'
index=${j}
set output 'm${j}.gif'


set border linewidth 0
unset tics
set lmargin screen 0.1
set rmargin screen 0.9
set tmargin screen 0.9
set bmargin screen 0.1


set palette model RGB defined ( \
                    0. 1 1 1,\
                    0.01 1 1 1,\
                    0.11 0 0 1,\
                    0.36 0.5 1 1,\
                    0.66 1 1 0,\
                    0.89 1 0 0,\
                    1 0.5 0 0 \
                )

		      
set pm3d map
do for [i=${start}:${filesCounter}]{echo "${files[${i}]}"; set title sprintf('T = %d dt', i); plot sprintf(${files[${i}]}) u 2:1:3 matrix with image}
EOD
