/*
 *
 */

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "seism_t.h"
#include "pthread.h"

int main(int argc, char** argv)
{
    int irank = -1;
    int nrank = -1;
    int provided = -1;
    seism_t_mpi_init(argc, argv, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);

    printf("RANK : %d, of %d : provided = %d\n", irank, nrank, provided);


    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}