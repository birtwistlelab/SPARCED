#!/bin/bash

echo hi > a.txt

# parse command-line arguments
if [[ $# != 1 ]]; then
	touch sweep.txt
  exit 1
fi

echo "$1" > sweep.txt

echo hello > b.txt
