/*
 *
 */

#include "seism_t.h"

/*
 *  initial vel_out pthread parameters
 */

int sgt_read_stations(int **tmpsta)
{

    FILE *fin;
    int i;
    printf("|    SGT: read stations...\n");
    fin = fopen(sgt_stationList,"r");
    if (fin == NULL)
    {
        printf("SGT: Cannot open station list\n");
        return -1;
    }
    fscanf(fin," %d", &sgt_numsta);
    printf("|    SGT: numsta = %d\n",sgt_numsta);
    *tmpsta = (int*)malloc(sizeof(int)*sgt_numsta*3);
    for(i=0;i<sgt_numsta;i++)
    {
        fscanf(fin,"%d %d %d", &(*tmpsta)[i*3], &(*tmpsta)[i*3+1], &(*tmpsta)[i*3+2]);
    }

    return 0;
}

/*
 *
 */

void sgt_out_init(struct sgt_input **sgt_args, int *sgt_opt, char *sgt_file,
                  int *areRegularReceivers, int *numsta, int *numsta_local,
                  int **sta, int **sta_indices, int *rank, int *coord,
                  int NZ, int nxt, int nyt, int nzt,
                  int *nt, int *NTISKP, int *WRITE_STEP,
                  int *cur_step_ptr, Grid3D lam, Grid3D mu,
                  pthread_mutex_t *ntiskp_mutex,
                  pthread_rwlock_t *veldata_rwlock,
                  pthread_cond_t *ntiskp_cond,
                  float **tmpbuf, MPI_File *fh)
{

}

/*
 *
 */

void sgt_out_exec(void *ptr)
{
    struct sgtout_input *sgtout_args;
    sgtout_args = (struct sgtout_input*)ptr;

    int *irank              = sgtout_args->irank;
    int *nt                 = sgtout_args->nt;
    int *NTISKP             = sgtout_args->NTISKP;
    int *WRITE_STEP         = sgtout_args->WRITE_STEP;
    int *cur_step_ptr       = sgtout_args->cur_step_ptr;
    float **tmpbuf          = sgtout_args->tmpbuf;
    MPI_File *fh            = sgtout_args->fh;
    MPI_Datatype filetype   = sgtout_args->filetype;
    MPI_Offset displacement = sgtout_args->displacement;

    int cur_step_local = *cur_step_ptr;
    int i, pos, ind, err;
    char filename[100];
    MPI_Status filestatus;
    char mpi_err_str[100];
    int mpi_err_str_l;

    for (i=0; i<sgt_numsta_local; i++)
    {
        pos = i*6;

        if (!i) printf("RANK: %d, sgt station ind = %d\n", *irank,i);
        if (!i) printf("\ttmpbuf(%d,%d,%d)=\n", sgt_sta[i*3], sgt_sta[i*3+1], sgt_sta[i*3+2]);
        if (!i) printf("\t%e,%e,%e,%e,%e,%e;\n", (*tmpbuf)[pos], (*tmpbuf)[pos+1], (*tmpbuf)[pos+2], (*tmpbuf)[pos+3], (*tmpbuf)[pos+4], (*tmpbuf)[pos+5]);

        sgt_buf[ind] = sgt_c1[i]*(*tmpbuf)[pos] + sgt_c2[i]*(*tmpbuf)[pos+1] + sgt_c2[i]*(*tmpbuf)[pos+2];
        ind += (*WRITE_STEP);

        sgt_buf[ind] = sgt_c2[i]*(*tmpbuf)[pos] + sgt_c1[i]*(*tmpbuf)[pos+1] + sgt_c2[i]*(*tmpbuf)[pos+2];
        ind += (*WRITE_STEP);

        sgt_buf[ind] = sgt_c2[i]*(*tmpbuf)[pos] + sgt_c2[i]*(*tmpbuf)[pos+1] + sgt_c1[i]*(*tmpbuf)[pos+2];
        ind += (*WRITE_STEP);

        sgt_buf[ind] = sgt_mu[i]*(*tmpbuf)[pos+3];
        ind += (*WRITE_STEP);

        sgt_buf[ind] = sgt_mu[i]*(*tmpbuf)[pos+4];
        ind += (*WRITE_STEP);

        sgt_buf[ind] = sgt_mu[i]*(*tmpbuf)[pos+5];
        ind += (*WRITE_STEP);
    }

    if ((cur_step_local/(*NTISKP))%(*WRITE_STEP) == 0)
    {
        pthread_mutex_lock(&mpi_mutex);

        if (*irank==0) printf("Writing sgt outputs...\n");
        sprintf(filename, "%s%07d", sgt_filenamebase, cur_step_local);
        err = MPI_File_open(sgtout_MCwriters,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,fh);
        if (err!=MPI_SUCCESS)
        {
            MPI_Error_string(err, mpi_err_str, &mpi_err_str_l);
            printf("%d) MPI_ERROR! %s\n",*irank,mpi_err_str);
        }
        err = MPI_File_set_view(*fh, displacement, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
        if (err!=MPI_SUCCESS)
        {
            MPI_Error_string(err, mpi_err_str, &mpi_err_str_l);
            printf("%d) MPI_ERROR! %s\n",*irank,mpi_err_str);
        }
        err = MPI_File_write_all(*fh, sgt_buf, sgt_numsta_local*6*(*WRITE_STEP), MPI_FLOAT, &filestatus);
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
        pthread_mutex_unlock(&mpi_mutex);
    }
}