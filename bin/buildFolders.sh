#!/bin/bash

numItems=$(ls -dq sweep* | wc -l)

for i in $( seq 0 $(( $numItems-1 )) )
    do
    rsync -avr --exclude="sweep*" --exclude="outputFolder*" "." "outputFolder$i"
    mv "sweep$i.txt" "outputFolder$i/"
done