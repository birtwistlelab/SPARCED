#!/bin/bash

numItems=$(ls -dq sweep* | wc -l)
numCopies=$1

for i in $( seq 0 $(( $numItems-1 )) )
    do
    for j in $( seq 1 ${numCopies} )
        do
        rsync -avr --exclude="sweep*" --exclude="outputFolder*" "." "outputFolder${i}copy${j}"
        cp "sweep$i.txt" "outputFolder${i}copy${j}"
    done
done