//
// Created by dmu on 8/9/17.
//

#ifndef SEISM_T_VEL_OUT_H
#define SEISM_T_VEL_OUT_H

#include "stdio.h"
#include "stdlib.h"
#include "mpi.h"

void seism_t_mpi_init(int argc, char** argv, int* provided);

void vel_out();

#endif //SEISM_T_VEL_OUT_H
