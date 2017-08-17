//
// Created by dmu on 8/9/17.
//

#ifndef SEISM_T_VEL_OUT_H
#define SEISM_T_VEL_OUT_H

#include "stdio.h"
#include "stdlib.h"
#include "memory.h"
#include "pthread.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "mpi.h"

pthread_t vel_out_thread;

pthread_attr_t   thread_attr;
pthread_mutex_t  ntiskp_mutex;
pthread_rwlock_t veldata_rwlock;
pthread_cond_t   ntiskp_cond;
pthread_mutex_t  vel_out_mutex;
pthread_cond_t   vel_out_cond;
pthread_mutex_t mpi_mutex;
pthread_cond_t mpi_cond;

void seism_t_mpi_init(int argc, char** argv, int* provided);

typedef float *restrict *restrict *restrict Grid3D;
typedef float *restrict Grid1D;
typedef int   *restrict PosInf;



pthread_t velout_thread;
MPI_Comm  velout_MCwriters;

struct velout_input
{
    int *irank;
    float *Bufx, *Bufy, *Bufz;
    float *tmpbuf;
    int *rec_nxt, *rec_nyt, *rec_nzt;
    int *nt, *NTISKP, *WRITE_STEP;
    int *cur_step_ptr;
    pthread_mutex_t *ntiskp_mutex;
    pthread_rwlock_t *veldata_rwlock;
    pthread_cond_t *ntiskp_cond;
    cudaStream_t *stream_IO;
    MPI_File *fh;
    MPI_Offset *displacement;
    MPI_Datatype *filetype;
    char *filename, *filenamebasex, *filenamebasey, *filenamebasez;
};

struct velout_input *velout_args;

void vel_out_param_init(struct velout_input **velout_args, int *irank,
                        float *Bufx, float *Bufy, float *Bufz, float *tmpbuf,
                        int *rec_nxt, int *rec_nyt, int *rec_nzt,
                        int *nt, int *NTISKP, int *WRITE_STEP,
                        int *cur_step_ntiskp,
                        pthread_mutex_t *ntiskp_mutex,
                        pthread_rwlock_t *veldata_rwlock,
                        pthread_cond_t *ntiskp_cond,
                        cudaStream_t *stream_IO, MPI_File *fh,
                        MPI_Offset *displacement, MPI_Datatype *filetype,
                        char *filename, char *filenamebasex, char *filenamebasey, char *filenamebasez);

void vel_out_pthread_func(void *ptr);

#endif //SEISM_T_VEL_OUT_H
