cmake_minimum_required(VERSION 3.0.2)

project(Nemoh Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "-cpp")
    # Other flags that may be useful: -r8 (force double-precision on all real numbers)
    add_definitions(-D_INTEL_FORTRAN_COMPILE_RULE_)
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "-cpp -DGNUFORT -ffree-line-length-none")
    # Other flags that may be useful: -fdefault-real-8 (force double-precision on all real numbers)
    add_definitions(-D_GNU_FORTRAN_COMPILE_RULE_)
endif()

set(COMMON_SRC_DIR ${PROJECT_SOURCE_DIR}/Common)
set(CORE_SRC_DIR ${PROJECT_SOURCE_DIR}/Solver/Core)

add_subdirectory(Mesh)
add_subdirectory(preProcessor)
add_subdirectory(postProcessor)
add_subdirectory(Solver)

add_subdirectory(QTF)

include(CTest)

add_subdirectory(TestCases)
