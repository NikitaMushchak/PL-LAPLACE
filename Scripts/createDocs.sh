#!/bin/bash

source ~/.bash_interface

printHeader "Generating docs..."

cd ./Sources/Sequential/

doxygen ../../Scripts/doxygen.html.2.0.txt

mv -v ./html ../../Docs