#!/bin/bash

# parse command-line arguments
if [[ $# != 1 ]]; then
	touch sweep.txt
  exit 1
fi

#enter paramsweep into sweep.txt
echo "$1" > sweep.txt
