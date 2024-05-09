# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_SPARCED.dir/SPARCEDPYTHON_wrap.cxx"
  "swig/SPARCED.py"
  )
endif()
