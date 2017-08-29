//
// Created by dmu on 8/9/17.
//

#ifndef SEISM_T_VEL_OUT_H
#define SEISM_T_VEL_OUT_H

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <pthread.h>
#include <mpi.h>
#include <unistd.h>
#include <math.h>

#define align 32
#define loop  1

pthread_t vel_thread;

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

MPI_Comm  velout_MCwriters;

struct velout_input
{
    int   *irank;
    int   *rec_nxt, *rec_nyt, *rec_nzt;
    int   *nt, *NTISKP, *WRITE_STEP;
    int   *cur_step;
    float *Bufx, *Bufy, *Bufz;
    char  *filenamebase, *filenamebasex, *filenamebasey, *filenamebasez;
    MPI_File     *fh;
    MPI_Offset   *displacement;
    MPI_Datatype *filetype;
};

struct velout_input *velout_args;

void vel_out_init(struct velout_input **velout_args, int *irank,
                  float *Bufx, float *Bufy, float *Bufz,
                  int *rec_nxt, int *rec_nyt, int *rec_nzt,
                  int *nt, int *NTISKP, int *WRITE_STEP,
                  int *cur_step, MPI_File *fh,
                  MPI_Offset *displacement, MPI_Datatype *filetype,
                  char *filenamebase, char *filenamebasex,
                  char *filenamebasey, char *filenamebasez);

void vel_out_exec(void *ptr);

void vel_out_pthread(struct velout_input* ptr);

//
// SGT FUNCTION
//

pthread_t sgt_thread;
int       sgt_numsta;
int       sgt_numsta_local;
float    *sgt_buf;
int      *sgt_sta;
char     *sgt_filenamebase;
char     *sgt_stationList;
int      *sgt_sta_indices;
float    *sgt_c1, *sgt_c2, *sgt_mu;
int       sgt_igreen;

struct sgt_input
{
    int *rank;
    int *nt;
    int *NTISKP;
    int *WRITE_STEP;
    int *cur_step_ptr;
    pthread_mutex_t *ntiskp_mutex;
    pthread_rwlock_t *veldata_rwlock;
    pthread_cond_t *ntiskp_cond;
    float **tmpbuf;
    MPI_File *fh;
    MPI_Offset displacement;
    MPI_Datatype filetype;
};



#endif //SEISM_T_VEL_OUT_H
