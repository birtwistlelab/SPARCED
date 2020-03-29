process start {
  '''
  pip -V > where.txt
  pip freeze > reqs.txt

  cd ..
  cd ..
  cd ..


  cd SPARCED
  pip install python-libsbml
  python3 SPARCED_ModelCreateWrite.py
  python3 RunModel.py
  '''
}
