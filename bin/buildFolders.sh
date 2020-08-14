#!/bin/bash

numItems=$(ls -dq sweep* | wc -l)

for i in $( seq 1 $numItems )
    do
    rsync -avr --exclude="sweep*" --exclude="outputFolder*" "." "outputFolder$i"
    cp "sweep$i.txt" "outputFolder$i"
done