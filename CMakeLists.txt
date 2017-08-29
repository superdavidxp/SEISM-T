cmake_minimum_required(VERSION 3.2)
project(SEISM_T)

set (CMAKE_C_STANDARD 99)

enable_language (Fortran)
find_package (Threads REQUIRED)
if (Threads_FOUND)
    message("")
    message("|    PTHREAD FOUND")
    message("|    PTHREAD = " ${CMAKE_THREAD_LIBS_INIT})
endif (Threads_FOUND)


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

#find_package(CUDA REQUIRED)
#if(CUDA_FOUND)
#    message("")
#    message("|    CUDA FOUND")
#    message("|    " ${CUDA_VERSION_STRING})
#    message("|    " ${CUDA_TOOLKIT_ROOT_DIR})
#    message("|    " ${CUDA_INCLUDE_DIRS})
#    message("|    " ${CUDA_LIBRARIES})
#endif(CUDA_FOUND)

SET(CMAKE_CXX_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
SET(CMAKE_C_FLAGS "-g -O0 -Wall -W -fprofile-arcs -ftest-coverage")
SET(CMAKE_EXE_LINKER_FLAGS "-fprofile-arcs -ftest-coverage")

set (CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

set (CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/seism_t)
set (CMAKE_LIBRARY_PATH   ${CMAKE_INSTALL_PREFIX}/lib)
set (CMAKE_INCLUDE_PATH   ${CMAKE_INSTALL_PREFIX}/include)

set (CMAKE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

message ("CMAKE_SOURCE_DIR     = " ${CMAKE_SOURCE_DIR})
message ("CMAKE_INSTALL_PREFIX = " ${CMAKE_INSTALL_PREFIX})
message ("CMAKE_LIBRARY_PATH   = " ${CMAKE_LIBRARY_PATH})
message (${CMAKE_INCLUDE_PATH})

include_directories(${CMAKE_SOURCE_DIR} ${MPI_C_INCLUDE_PATH})

add_library (seism_t SHARED src/vel_out.c src/seism_t.h src/common.c)
target_link_libraries(seism_t ${CMAKE_THREAD_LIBS_INIT} ${MPI_C_LIBRARIES})

INSTALL (TARGETS seism_t DESTINATION ${CMAKE_LIBRARY_PATH})
INSTALL (DIRECTORY src/ DESTINATION ${CMAKE_INCLUDE_PATH}
         FILES_MATCHING PATTERN "*.h")

add_executable(examples_01 examples/examples_01.c)
target_link_libraries(examples_01 seism_t ${CMAKE_THREAD_LIBS_INIT} ${MPI_C_LIBRARIES} -lm)

add_executable(examples_02 examples/examples_02.c)
target_link_libraries(examples_02 seism_t ${CMAKE_THREAD_LIBS_INIT} ${MPI_C_LIBRARIES} -lm)

#enable_testing()
#add_test(EXAMPLE_01 examples)

