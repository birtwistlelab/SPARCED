cmake_minimum_required(VERSION 2.8)

if(POLICY CMP0060)
  cmake_policy(SET CMP0060 NEW)
endif(POLICY CMP0060)
if(POLICY CMP0065)
  cmake_policy(SET CMP0065 NEW)
endif(POLICY CMP0065)

project(SPARCED_Brep)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CheckCXXCompilerFlag)
set(MY_CXX_FLAGS -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable)
foreach(FLAG ${MY_CXX_FLAGS})
    unset(CUR_FLAG_SUPPORTED CACHE)
    CHECK_CXX_COMPILER_FLAG(${FLAG} CUR_FLAG_SUPPORTED)
    if(${CUR_FLAG_SUPPORTED})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
    endif()
endforeach(FLAG)

find_package(Amici 0.10.13 HINTS ${CMAKE_CURRENT_LIST_DIR}/../../build)
message(STATUS "Found AMICI ${Amici_DIR}")

set(MODEL_DIR ${CMAKE_CURRENT_LIST_DIR})

set(SRC_LIST_LIB SPARCED_Brep_J.cpp 
SPARCED_Brep_JB.cpp 
SPARCED_Brep_JDiag.cpp 
SPARCED_Brep_JSparse.cpp 
SPARCED_Brep_JSparseB.cpp 
SPARCED_Brep_Jy.cpp 
SPARCED_Brep_dJydsigmay.cpp 
SPARCED_Brep_dJydy.cpp 
SPARCED_Brep_dwdp.cpp 
SPARCED_Brep_dwdx.cpp 
SPARCED_Brep_dxdotdw.cpp 
SPARCED_Brep_dxdotdp.cpp 
SPARCED_Brep_dydx.cpp 
SPARCED_Brep_dydp.cpp 
SPARCED_Brep_dsigmaydp.cpp 
SPARCED_Brep_sigmay.cpp 
SPARCED_Brep_w.cpp 
SPARCED_Brep_x0.cpp 
SPARCED_Brep_x0_fixedParameters.cpp 
SPARCED_Brep_sx0.cpp 
SPARCED_Brep_sx0_fixedParameters.cpp 
SPARCED_Brep_xdot.cpp 
SPARCED_Brep_y.cpp 
SPARCED_Brep_x_rdata.cpp 
SPARCED_Brep_total_cl.cpp 
SPARCED_Brep_x_solver.cpp 
SPARCED_Brep_JSparse_colptrs.cpp
SPARCED_Brep_JSparse_rowvals.cpp 
SPARCED_Brep_JSparseB_colptrs.cpp
SPARCED_Brep_JSparseB_rowvals.cpp 
SPARCED_Brep_dJydy_colptrs.cpp
SPARCED_Brep_dJydy_rowvals.cpp 
SPARCED_Brep_dwdx_colptrs.cpp
SPARCED_Brep_dwdx_rowvals.cpp 
SPARCED_Brep_dxdotdw_colptrs.cpp
SPARCED_Brep_dxdotdw_rowvals.cpp 
${MODEL_DIR}/wrapfunctions.cpp
)

add_library(${PROJECT_NAME} ${SRC_LIST_LIB})
add_library(model ALIAS ${PROJECT_NAME})

# This option can be helpful when using the Intel compiler and compilation of
# wrapfunctions.cpp fails due to insufficient memory.
option(ENABLE_WRAPFUNCTIONS_OPTIMIZATIONS "Enable compiler optimizations for wrapfunctions.cpp?" ON)
if(NOT ENABLE_WRAPFUNCTIONS_OPTIMIZATIONS)
    set_source_files_properties(wrapfunctions.cpp PROPERTIES COMPILE_FLAGS -O0)
endif()

target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(${PROJECT_NAME}
    PUBLIC Upstream::amici
)

set(SRC_LIST_EXE main.cpp)

add_executable(simulate_${PROJECT_NAME} ${SRC_LIST_EXE})

target_link_libraries(simulate_${PROJECT_NAME} ${PROJECT_NAME})

if($ENV{ENABLE_GCOV_COVERAGE})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()

## SWIG
option(ENABLE_SWIG "Build swig/python library?" ON)
if(ENABLE_SWIG)
    if(NOT(${CMAKE_VERSION} VERSION_LESS 3.8))
        add_subdirectory(swig)
    else()
        message(WARNING "Unable to build SWIG interface, upgrade CMake to >=3.8.")
    endif()
endif()


# <Export cmake configuration>
include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}Targets
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
export(EXPORT ${PROJECT_NAME}Targets FILE ${PROJECT_NAME}Config.cmake
    NAMESPACE Upstream::
    )
# </Export cmake configuration>

