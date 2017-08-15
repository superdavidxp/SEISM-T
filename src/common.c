/*
 *
 */

#include "seism_t.h"

void seism_t_mpi_init(int argc, char** argv, int* provided)
{
    int err;
    err=MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, provided);
    if (err != 0)
    {
        printf("|    ERROR: seism_t_mpi_init failed.\n");
        exit(1);
    }
    else
    {
        printf("|    seism_t_mpi_init succesful. provided = %d\n", *provided);
    }
}

