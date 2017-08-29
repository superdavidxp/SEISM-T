/*
 *
 */

#include "seism_t.h"

/*
 *  initial vel_out pthread parameters
 */

void vel_out_init(struct velout_input **velout_args, int *irank,
                  float *Bufx, float *Bufy, float *Bufz,
                  int *rec_nxt, int *rec_nyt, int *rec_nzt,
                  int *nt, int *NTISKP, int *WRITE_STEP,
                  int *cur_step, MPI_File *fh,
                  MPI_Offset *displacement, MPI_Datatype *filetype,
                  char *filenamebase, char *filenamebasex,
                  char *filenamebasey, char *filenamebasez)
{
    int err = -1;

    *velout_args = (struct velout_input*)malloc(sizeof(struct velout_input));

    (*velout_args)->irank = irank;
    (*velout_args)->Bufx  = Bufx;
    (*velout_args)->Bufy  = Bufy;
    (*velout_args)->Bufz  = Bufz;
    (*velout_args)->rec_nxt = rec_nxt;
    (*velout_args)->rec_nyt = rec_nyt;
    (*velout_args)->rec_nzt = rec_nzt;
    (*velout_args)->nt      = nt;
    (*velout_args)->NTISKP  = NTISKP;
    (*velout_args)->WRITE_STEP  = WRITE_STEP;
    (*velout_args)->cur_step = cur_step;
    (*velout_args)->fh  = fh;
    (*velout_args)->displacement  = displacement;
    (*velout_args)->filetype  = filetype;
    (*velout_args)->filenamebase  = filenamebase;
    (*velout_args)->filenamebasex = filenamebasex;
    (*velout_args)->filenamebasey = filenamebasey;
    (*velout_args)->filenamebasez = filenamebasez;
}

void vel_out_pthread(struct velout_input* ptr)
{
    int err    = -1;
    int *irank = ptr->irank;

    err = pthread_create(&vel_thread, &thread_attr, (void*) &vel_out_exec, (void*) ptr);
    if (err)
    {
        printf("ERROR: RANK = %d: Cannot create velout thread! Error Code = %d\n", *irank, err);
        exit(1);
    }
}

void vel_out_exec(void *ptr)
{
    struct velout_input* velout_args;
    velout_args     = (struct velout_input*)ptr;
    int *irank      = velout_args->irank;
    int *rec_nxt    = velout_args->rec_nxt;
    int *rec_nyt    = velout_args->rec_nyt;
    int *rec_nzt    = velout_args->rec_nzt;
    int *nt         = velout_args->nt;
    int *NTISKP     = velout_args->NTISKP;
    int *WRITE_STEP = velout_args->WRITE_STEP;
    int *cur_step   = velout_args->cur_step;
    float *Bufx     = velout_args->Bufx;
    float *Bufy     = velout_args->Bufy;
    float *Bufz     = velout_args->Bufz;

    MPI_File     *fh = velout_args->fh;
    MPI_Offset   *displacement = velout_args->displacement;
    MPI_Datatype *filetype = velout_args->filetype;
    MPI_Status   filestatus;

    char *filename  = velout_args->filenamebase;
    char *filenamebasex = velout_args->filenamebasex;
    char *filenamebasey = velout_args->filenamebasey;
    char *filenamebasez = velout_args->filenamebasez;

    long int idtmp, num_bytes;
    int err;
    int cur_step_local;
    cur_step_local = *cur_step;

    pthread_mutex_lock( &vel_out_mutex );

    if ( fmod(((*cur_step)/(*NTISKP)), *WRITE_STEP )==0 )
    {
        if (*irank == 0) printf("|    RANK %d: In vel_out_exec! cur_step=%d\n",*irank, cur_step_local);

        //
        // Build a pthd buffer

        sprintf(filename, "%s%07d", filenamebasex, cur_step_local);
        err = MPI_File_open(velout_MCwriters, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, fh);
        err = MPI_File_set_view(*fh, *displacement, MPI_FLOAT, *filetype, "native", MPI_INFO_NULL);
        err = MPI_File_write_all(*fh, Bufx, (*rec_nxt) * (*rec_nyt) * (*rec_nzt) * (*WRITE_STEP), MPI_FLOAT,
                                 &filestatus);
        err = MPI_File_close(fh);

        sprintf(filename, "%s%07d", filenamebasey, cur_step_local);
        err = MPI_File_open(velout_MCwriters, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, fh);
        err = MPI_File_set_view(*fh, *displacement, MPI_FLOAT, *filetype, "native", MPI_INFO_NULL);
        err = MPI_File_write_all(*fh, Bufy, (*rec_nxt) * (*rec_nyt) * (*rec_nzt) * (*WRITE_STEP), MPI_FLOAT,
                                 &filestatus);
        err = MPI_File_close(fh);

        sprintf(filename, "%s%07d", filenamebasez, cur_step_local);
        err = MPI_File_open(velout_MCwriters, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, fh);
        err = MPI_File_set_view(*fh, *displacement, MPI_FLOAT, *filetype, "native", MPI_INFO_NULL);
        err = MPI_File_write_all(*fh, Bufz, (*rec_nxt) * (*rec_nyt) * (*rec_nzt) * (*WRITE_STEP), MPI_FLOAT,
                                 &filestatus);
        err = MPI_File_close(fh);

        // add one file sync to make sure the file is written to the disk before program complete

        if (*irank == 0) printf("|    Outputs are written successfully! sur_step = %d\n", *cur_step);

    } // if ( fmod((cur_step/NTISKP), WRITE_STEP )==0 )

    pthread_mutex_unlock( &vel_out_mutex );

    pthread_exit(0);

    return;
}

void vel_out_finalize()
{
    pthread_mutex_lock( &vel_out_mutex );

    pthread_join(vel_thread, NULL);
    printf("|    vel_thread ended.\n");

    pthread_mutex_unlock( &vel_out_mutex );

    return;
}

//char mpi_err_str[100];
//int mpi_err_str_l;
//    if (err!=MPI_SUCCESS)
//    {
//        MPI_Error_string(err, mpi_err_str, &mpi_err_str_l);
//        printf("%d) MPI_ERROR! %s\n",*irank,mpi_err_str);
//    }