process start {
  '''
  pip3 freeze > reqs.txt

  cd ..
  cd ..
  cd ..

  sleep 10000
  cd SPARCED
  python3 SPARCED_ModelCreateWrite.py
  python3 RunModel.py
  '''
}
