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

    pthread_rwlock_rdlock(veldata_rwlock);

    for (i=0; i<sgt_numsta_local; i++)
    {
        pos = i*6;

        if (!i) printf("%d) sgt station ind=%d\n", *rank,i);
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

    pthread_rwlock_unlock(veldata_rwlock);

    if ((cur_step_local/(*NTISKP))%(*WRITE_STEP) == 0)
    {
        pthread_mutex_lock(&mpi_mutex);

        if (*rank==0) printf("Writing sgt outputs...\n");
        sprintf(filename, "%s%07d", sgt_filenamebase, cur_step_local);
        err = MPI_File_open(sgt_MCW,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,fh);
        if (err!=MPI_SUCCESS)
        {
            MPI_Error_string(err, mpi_err_str, &mpi_err_str_l);
            printf("%d) MPI_ERROR! %s\n",*rank,mpi_err_str);
        }
        err = MPI_File_set_view(*fh, displacement, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
        if (err!=MPI_SUCCESS)
        {
            MPI_Error_string(err, mpi_err_str, &mpi_err_str_l);
            printf("%d) MPI_ERROR! %s\n",*rank,mpi_err_str);
        }
        err = MPI_File_write_all(*fh, sgt_buf, sgt_numsta_local*6*(*WRITE_STEP), MPI_FLOAT, &filestatus);
        if (err!=MPI_SUCCESS)
        {
            MPI_Error_string(err, mpi_err_str, &mpi_err_str_l);
            printf("%d) MPI_ERROR! %s\n",*rank,mpi_err_str);
        }
        err = MPI_File_close(fh);
        if (err!=MPI_SUCCESS)
        {
            MPI_Error_string(err, mpi_err_str, &mpi_err_str_l);
            printf("%d) MPI_ERROR! %s\n",*rank,mpi_err_str);
        }
        pthread_mutex_unlock(&mpi_mutex);
    }
}