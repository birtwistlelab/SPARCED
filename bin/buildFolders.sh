#!/bin/bash

numItems=$(ls -dq sweep* | wc -l)

for i in $( seq 0 $(( $numItems-1 )) )
    do
    rsync -avr --exclude="sweep*" --exclude="outputFolder*" "." "outputFolder$i"
    cp "sweep$i.txt" "outputFolder$i"
done