/*
 *
 */

#include "seism_t.h"

/*
 *  initial vel_out pthread parameters
 */
void vel_out_init()
{
    int err = -1;
    err=pthread_mutex_init(&vel_out_mutex,NULL);
    if (err != 0)
    {
        printf("|    ERROR: pthread_mutex_init failed.\n");
        exit(1);
    }
    else
    {
        printf("|    pthread_mutex_init succesful.\n");
    }
    err=pthread_cond_init(&vel_out_cond,NULL);
    if (err != 0)
    {
        printf("|    ERROR: pthread_cond_init failed.\n");
        exit(1);
    }
    else
    {
        printf("|    pthread_cond_init succesful.\n");
    }

}

/*
 *
 */

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
                        char *filename, char *filenamebasex, char *filenamebasey, char *filenamebasez)
{
    int err = -1;

    *velout_args = (struct velout_input*)malloc(sizeof(struct velout_input));

    (*velout_args)->irank = irank;
    (*velout_args)->Bufx  = Bufx;
    (*velout_args)->Bufy  = Bufy;
    (*velout_args)->Bufz  = Bufz;
    (*velout_args)->tmpbuf  = tmpbuf;
    (*velout_args)->rec_nxt = rec_nxt;
    (*velout_args)->rec_nyt = rec_nyt;
    (*velout_args)->rec_nzt = rec_nzt;
    (*velout_args)->nt      = nt;
    (*velout_args)->NTISKP  = NTISKP;
    (*velout_args)->WRITE_STEP  = WRITE_STEP;
    (*velout_args)->cur_step_ptr = cur_step_ntiskp;
    (*velout_args)->ntiskp_mutex  = ntiskp_mutex;
    (*velout_args)->veldata_rwlock  = veldata_rwlock;
    (*velout_args)->ntiskp_cond  = ntiskp_cond;
    (*velout_args)->stream_IO   = stream_IO;
    (*velout_args)->fh  = fh;
    (*velout_args)->displacement  = displacement;
    (*velout_args)->filetype  = filetype;
    (*velout_args)->filename  = filename;
    (*velout_args)->filenamebasex = filenamebasex;
    (*velout_args)->filenamebasey = filenamebasey;
    (*velout_args)->filenamebasez = filenamebasez;
}

void vel_out_pthread_func(void *ptr)
{
    struct velout_input* velout_args;
    velout_args       = (struct velout_input*)ptr;
    int *irank        = velout_args->irank;
    float *Bufx       = velout_args->Bufx;
    float *Bufy       = velout_args->Bufy;
    float *Bufz       = velout_args->Bufz;
    float *tmpbuf     = velout_args->tmpbuf;
    int *rec_nxt      = velout_args->rec_nxt;
    int *rec_nyt      = velout_args->rec_nyt;
    int *rec_nzt      = velout_args->rec_nzt;
    int *nt           = velout_args->nt;
    int *NTISKP       = velout_args->NTISKP;
    int *WRITE_STEP   = velout_args->WRITE_STEP;
    int *cur_step_ptr = velout_args->cur_step_ptr;

    pthread_mutex_t *ntiskp_mutex = velout_args->ntiskp_mutex;
    pthread_rwlock_t *veldata_rwlock = velout_args->veldata_rwlock;
    pthread_cond_t *ntiskp_cond = velout_args->ntiskp_cond;

    cudaStream_t *stream_IO  = velout_args->stream_IO;
    MPI_File *fh = velout_args->fh;

    MPI_Offset *displacement = velout_args->displacement;
    MPI_Datatype *filetype = velout_args->filetype;
    MPI_Status filestatus;

    char *filename  = velout_args->filename;
    char *filenamebasex = velout_args->filenamebasex;
    char *filenamebasey = velout_args->filenamebasey;
    char *filenamebasez = velout_args->filenamebasez;

    long int idtmp, num_bytes;
    int err;
    int cur_step_local;
    cudaError_t cerr;
    cur_step_local = *cur_step_ptr;
    char mpi_err_str[100];
    int mpi_err_str_l;

    printf("%d) In velout! cur_step=%d\n",*irank,cur_step_local);

    if(*irank==0) printf("velout cur_step=%d\n",cur_step_local);
    if(*irank==0) printf("Args: NTISKP,WRITE_STEP=%d,%d\n",*NTISKP,*WRITE_STEP);

    while (cur_step_local+(*NTISKP)<=(*nt))
    {
        if(*irank==0)
        printf("velout waiting for ntiskp_mutex inside loop cur_step=%d\n", cur_step_local);
        pthread_mutex_lock(ntiskp_mutex);
        while(0!=pthread_cond_wait(ntiskp_cond, ntiskp_mutex));
        cur_step_local = *cur_step_ptr;
        pthread_mutex_unlock(ntiskp_mutex);
//        if(*irank==0)
//        printf("velout inside loop cur_step=%d\n",cur_step_local);
//
//        num_bytes = sizeof(float)*(*rec_nxt)*(*rec_nyt)*(*rec_nzt);
//        idtmp = ((cur_step_local/(*NTISKP)+(*WRITE_STEP)-1)%(*WRITE_STEP));
//        idtmp = idtmp*(*rec_nxt)*(*rec_nyt)*(*rec_nzt);
//
//        cudaStreamSynchronize(*stream_IO);
//
//        cerr = cudaGetLastError();
//        if (cerr != cudaSuccess) printf("CUDA KERNEL ERROR! velout after stream_IO is synchronized:%s\n",cudaGetErrorString(cerr));

        pthread_rwlock_rdlock(veldata_rwlock);

//        memcpy(Bufx+idtmp, tmpbuf, num_bytes);
//        memcpy(Bufy+idtmp, tmpbuf+num_bytes/sizeof(float), num_bytes);
//        memcpy(Bufz+idtmp, tmpbuf+num_bytes/sizeof(float)*2, num_bytes);

        pthread_rwlock_unlock(veldata_rwlock);

        if((cur_step_local/(*NTISKP))%(*WRITE_STEP) == 0)
        {
            pthread_mutex_lock(&mpi_mutex);

            if (*irank==0) printf("Writing outputs...\n");
            sprintf(filename, "%s%07d", filenamebasex, cur_step_local);
            err = MPI_File_open(velout_MCwriters,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,fh);
            if (err!=MPI_SUCCESS)
            {
                MPI_Error_string(err, mpi_err_str, &mpi_err_str_l);
                printf("%d) MPI_ERROR! %s\n",*irank,mpi_err_str);
            }
            err = MPI_File_set_view(*fh, *displacement, MPI_FLOAT, *filetype, "native", MPI_INFO_NULL);
            if (err!=MPI_SUCCESS)
            {
                MPI_Error_string(err, mpi_err_str, &mpi_err_str_l);
                printf("%d) MPI_ERROR! %s\n",*irank,mpi_err_str);
            }
            err = MPI_File_write_all(*fh, Bufx, (*rec_nxt)*(*rec_nyt)*(*rec_nzt)*(*WRITE_STEP), MPI_FLOAT, &filestatus);
            if (err!=MPI_SUCCESS)
            {
                MPI_Error_string(err, mpi_err_str, &mpi_err_str_l);
                printf("%d) MPI_ERROR! %s\n",*irank,mpi_err_str);
            }
            err = MPI_File_close(fh);
            if (err!=MPI_SUCCESS)
            {
                MPI_Error_string(err, mpi_err_str, &mpi_err_str_l);
                printf("%d) MPI_ERROR! %s\n",*irank,mpi_err_str);
            }
            sprintf(filename, "%s%07d", filenamebasey, cur_step_local);
            err = MPI_File_open(velout_MCwriters,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,fh);
            err = MPI_File_set_view(*fh, *displacement, MPI_FLOAT, *filetype, "native", MPI_INFO_NULL);
            err = MPI_File_write_all(*fh, Bufx, (*rec_nxt)*(*rec_nyt)*(*rec_nzt)*(*WRITE_STEP), MPI_FLOAT, &filestatus);
            err = MPI_File_close(fh);
            sprintf(filename, "%s%07d", filenamebasez, cur_step_local);
            err = MPI_File_open(velout_MCwriters,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,fh);
            err = MPI_File_set_view(*fh, *displacement, MPI_FLOAT, *filetype, "native", MPI_INFO_NULL);
            err = MPI_File_write_all(*fh, Bufx, (*rec_nxt)*(*rec_nyt)*(*rec_nzt)*(*WRITE_STEP), MPI_FLOAT, &filestatus);
            err = MPI_File_close(fh);
            if (*irank==0) printf("Outputs are written successfully!\n");

            pthread_mutex_unlock(&mpi_mutex);

        } // if((cur_step_local/(*NTISKP))%(*WRITE_STEP) == 0)

    } // while (cur_step_local+(*NTISKP)<=(*nt))

    pthread_exit(0);

    return;
}