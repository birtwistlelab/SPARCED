#!/bin/bash

# parse command-line arguments
if [[ $# != 1 ]]; then
	touch sweep.txt
else
	#enter paramsweep into sweep.txt
	echo "$1" > sweep.txt
fi
