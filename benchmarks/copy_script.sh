#!/bin/bash

# Get the absolute path of the script directory
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Find all subdirectories in the current directory
SUB_DIRS=$(find . -mindepth 1 -maxdepth 1 -type d -path ./src -prune)

# Copy the file into each subdirectory
for dir in $SUB_DIRS; do
    cp "$SCRIPT_DIR/src/runModel_unitTest.py" "$dir/scripts/"
done

echo "Files copied successfully into subdirectories."
