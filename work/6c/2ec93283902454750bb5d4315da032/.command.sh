#!/bin/bash -e
pip3 freeze > reqs.txt

cd ..
cd ..
cd ..

cd SPARCED
python3 SPARCED_ModelCreateWrite.py
python3 RunModel.py
