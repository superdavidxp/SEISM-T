/*
 *
 */

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "seism_t.h"
#include "pthread.h"

void *print_message_function( void *ptr )
{
    char *message;
    message = (char *) ptr;
    printf("%s \n", message);
}

int main(int argc, char** argv)
{
    int irank, nrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);

    printf("|    RANK : %d, TEST STARTS\n", irank);

    pthread_t* thread =  (pthread_t*) malloc (2 * sizeof(pthread_t));

    char *message1 = "Thread main";
    char *message2 = "Thread vel_out";
    int  iret1, iret2;

    /* Create independent threads each of which will execute function */

    iret1 = pthread_create( &thread[0], NULL, print_message_function, (void*) message1);
    iret2 = pthread_create( &thread[1], NULL, print_message_function, (void*) message2);

    /* Wait till threads are complete before main continues. Unless we  */
    /* wait we run the risk of executing an exit which will terminate   */
    /* the process and all threads before the threads have completed.   */

    pthread_join( thread[0], NULL);
    pthread_join( thread[1], NULL);

    printf("Thread main returns: %d\n",iret1);
    printf("Thread vel_out returns: %d\n",iret2);

    MPI_Barrier(MPI_COMM_WORLD);

    return 0;
}