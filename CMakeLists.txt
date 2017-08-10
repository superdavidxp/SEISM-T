cmake_minimum_required(VERSION 3.2)
project(SEISM_T)

set (CMAKE_C_STANDARD 99)

enable_language (Fortran)

find_package (MPI REQUIRED)
if (MPI_FOUND)
    message("")
    message("|    MPI FOUND")
    message("|    MPIEXEC                  = " ${MPIEXEC})
    message("|    MPI_C_COMPILER           = " ${MPI_C_COMPILER})
    message("|    MPI_C_INCLUDE_PATH       = " ${MPI_C_INCLUDE_PATH})
    message("|    MPI_C_LIBRARIES          = " ${MPI_C_LIBRARIES})
    message("|    MPI_CXX_COMPILER         = " ${MPI_CXX_COMPILER})
    message("|    MPI_CXX_INCLUDE_PATH     = " ${MPI_CXX_INCLUDE_PATH})
    message("|    MPI_CXX_LIBRARIES        = " ${MPI_CXX_LIBRARIES})
    message("|    MPI_Fortran_COMPILER     = " ${MPI_Fortran_COMPILER})
    message("|    MPI_Fortran_INCLUDE_PATH = " ${MPI_Fortran_INCLUDE_PATH})
    message("|    MPI_Fortran_LIBRARIES    = " ${MPI_Fortran_LIBRARIES})
    message("|    MPI_INCLUDE_PATH         = " ${MPI_INCLUDE_PATH})
    message("|    MPI_LIBRARY              = " ${MPI_LIBRARY})
endif (MPI_FOUND)

set (CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

set (CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/seism_t)
set (CMAKE_LIBRARY_PATH   ${CMAKE_INSTALL_PREFIX}/lib)
set (CMAKE_INCLUDE_PATH   ${CMAKE_INSTALL_PREFIX}/include)

set (CMAKE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

message ("CMAKE_SOURCE_DIR     = " ${CMAKE_SOURCE_DIR})
message ("CMAKE_INSTALL_PREFIX = " ${CMAKE_INSTALL_PREFIX})
message ("CMAKE_LIBRARY_PATH   = " ${CMAKE_LIBRARY_PATH})
message (${CMAKE_INCLUDE_PATH})

add_library (seism_t SHARED src/vel_out.c src/seism_t.h)

INSTALL (TARGETS seism_t DESTINATION ${CMAKE_LIBRARY_PATH})
INSTALL (DIRECTORY src/ DESTINATION ${CMAKE_INCLUDE_PATH}
         FILES_MATCHING PATTERN "*.h")

include_directories(${CMAKE_SOURCE_DIR})
add_executable(example_01 examples/examples.c)
target_link_libraries(example_01 seism_t)

enable_testing()
add_test(EXAMPLE_01 example_01)

