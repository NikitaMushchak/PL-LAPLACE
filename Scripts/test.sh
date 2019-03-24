#!/bin/bash

CXX='g++'
CXFLAGS='-O3 -std=c++11 -I../Include/catch2 -I../Sources/Sequential'

mainName=tests-main
mainLink=./${mainName}.o
testOut=./tests.out
testFiles="tests-rhs3D.cc"

cd ./Tests

rm -f ${testOut}

if ! [ -f ${mainLink} ]; then
    echo "Compiling main..."
    ${CXX} ${CXFLAGS} -o ${mainLink} ${mainName}.cc -c
fi

echo "Compiling tests..."
${CXX} ${CXFLAGS} -o ${testOut} ${mainLink} ${testFiles} && ${testOut} -r compact