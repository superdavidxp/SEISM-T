/*
 *
 */

#include "seism_t.h"

/*
 *  initial vel_out pthread parameters
 */

int sgt_readStations(int **tmpsta)
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

void sgt_initialize(struct sgt_input **sgt_args, int *sgt_opt, char *sgt_file,
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
    FILE *fin;
    char opt[50];
    int i, ind;
    int *t_sgt_sta;
    int *t_sgt_sta_indices;
    int nbgx, nbgy, nbgz, nedx, nedy, nedz;
    int ind0, ind1, ind2;
    float tmplam, tmpmu;
    MPI_Offset displacement;
    int *ones;
    MPI_Aint *dispArray;
    MPI_Datatype filetype;

    *sgt_opt = 1;
    nbgx = coord[0]*nxt + 1;
    nbgy = coord[1]*nyt + 1;
    nbgz = 1;
    nedx = nbgx + nxt - 1;
    nedy = nbgy + nyt - 1;
    nedz = nbgz + nzt - 1;

    if (*rank == 0)
    {
        i = strcmp(sgt_file,"");
        if (i==0)
        {
            printf("SGT: no sgt parameters file is given!\n");
            *sgt_opt = 0;
        }
        else
        {
            fin = fopen(sgt_file,"r");
            if (fin==NULL)
            {
                printf("SGT: cannot open sgt-file!\n");
                *sgt_opt = 0;
            }
            else
            {
                sgt_filenamebase = (char*)malloc(sizeof(char)*50);
                while(fscanf(fin," %s",opt)!=EOF)
                {
                    if(strcmp(opt,"station_list")==0)
                    {
                        sgt_stationList = (char*)malloc(sizeof(char)*50);
                        fscanf(fin," %s",sgt_stationList);
                        printf("SGT OPTION: station_list: %s\n",sgt_stationList);
                    }
                    else if(strcmp(opt,"igreen")==0)
                    {
                        fscanf(fin," %d",&sgt_igreen);
                        printf("SGT OPTION: igreen: %d\n",sgt_igreen);
                    }
                    else if(strcmp(opt,"filenamebase")==0)
                    {
                        sgt_filenamebase = (char*)malloc(sizeof(char)*50);
                        fscanf(fin," %s",sgt_filenamebase);
                        printf("SGT OPTION: filenamebase: %s\n",sgt_filenamebase);
                    }
                    else{
                        printf("SGT: Unknown parameter: %s\n",opt);
                        free(sgt_stationList);
                        *sgt_opt = 0;
                    }
                }
                fclose(fin);
            }
        }
        if (*sgt_opt)
        {
            if(sgt_readStations(&t_sgt_sta))
            {
                *sgt_opt = 0;
            }
        }
    }
    printf("%d) SGT: share parameters\n",*rank);
    MPI_Comm_dup(MPI_COMM_WORLD, &sgt_MCW);
    pthread_mutex_lock(&mpi_mutex);
    printf("%d) SGT: bcast parameters\n",*rank);
    MPI_Bcast(sgt_opt,  1, MPI_INT,  0, sgt_MCW);
    pthread_mutex_unlock(&mpi_mutex);
    printf("%d) SGT: sgt_opt received: %d\n",*rank,*sgt_opt);
    if (*sgt_opt == 0)
    {
        return;
    }
    pthread_mutex_lock(&mpi_mutex);
    MPI_Bcast(&sgt_numsta, 1, MPI_INT, 0, sgt_MCW);
    printf("%d) SGT: sgt_numsta received: %d\n",*rank,sgt_numsta);
    if(*rank != 0) t_sgt_sta = (int*)malloc(sizeof(int)*sgt_numsta*3);
    MPI_Bcast(t_sgt_sta, sgt_numsta*3, MPI_INT, 0, sgt_MCW);
    printf("%d) SGT: sgt_sta received: (%d,%d,%d),...\n",*rank,
           t_sgt_sta[0], t_sgt_sta[1], t_sgt_sta[2]);
    if(*rank != 0) sgt_filenamebase = (char*)malloc(sizeof(char)*50);
    MPI_Bcast(sgt_filenamebase, 50, MPI_CHAR, 0, sgt_MCW);
    printf("%d) SGT: sgt_filenamebase received: %s\n",*rank,sgt_filenamebase);
    pthread_mutex_unlock(&mpi_mutex);
    printf("%d) SGT parameters are received successfully\n",*rank);
    sgt_numsta_local = 0;
    t_sgt_sta_indices = (int*)malloc(sizeof(int)*sgt_numsta);
    for(i=0;i<sgt_numsta;i++){
        t_sgt_sta[i*3+2] = NZ+1-t_sgt_sta[i*3+2];
        if(t_sgt_sta[i*3] >= nbgx && t_sgt_sta[i*3] <= nedx
           && t_sgt_sta[i*3+1] >= nbgy && t_sgt_sta[i*3+1] <= nedy
           && t_sgt_sta[i*3+2] >= nbgz && t_sgt_sta[i*3+2] <= nedz){
            t_sgt_sta_indices[sgt_numsta_local] = i;
            sgt_numsta_local++;
        }
    }
    if(sgt_numsta_local){
        printf("SGT: rank=%d, numsta=%d\n",*rank,sgt_numsta_local);
        sgt_sta = (int*)malloc(sizeof(int)*sgt_numsta_local*3);
        *sta = (int*)malloc(sizeof(int)*sgt_numsta_local*3);
        sgt_sta_indices = (int*)malloc(sizeof(int)*sgt_numsta_local);
        *sta_indices = (int*)malloc(sizeof(int)*sgt_numsta_local);

        MPI_Type_contiguous((*WRITE_STEP)*6, MPI_FLOAT, &filetype);
        MPI_Type_commit(&filetype);
        ones = (int*)malloc(sizeof(int)*sgt_numsta_local);
        dispArray = (MPI_Aint*)malloc(sizeof(MPI_Aint)*sgt_numsta_local);
        displacement = 0;

        for(i=0;i<sgt_numsta_local;i++){
            ind = t_sgt_sta_indices[i];
            sgt_sta_indices[i] = ind;
            sgt_sta[i*3] = t_sgt_sta[ind*3] - nbgx + 1;
            sgt_sta[i*3+1] = t_sgt_sta[ind*3+1] - nbgy + 1;
            sgt_sta[i*3+2] = t_sgt_sta[ind*3+2] - nbgz + 1;
            // sta and sta_indices may contain (in the future) stations
            //  for other supported pthread modules
            (*sta)[i*3]   = sgt_sta[i*3];
            (*sta)[i*3+1] = sgt_sta[i*3+1];
            (*sta)[i*3+2] = sgt_sta[i*3+2];
            (*sta_indices)[i] = sgt_sta_indices[i];

            dispArray[i] = sizeof(float)*6*(*WRITE_STEP)*ind - displacement;
            if(!i){ displacement = dispArray[0]; dispArray[0] = 0;}
            ones[i] = 1;
        }
    }
    free(t_sgt_sta);
    free(t_sgt_sta_indices);

    MPI_Type_create_hindexed(sgt_numsta_local, ones, dispArray, filetype, &filetype);
    MPI_Type_commit(&filetype);
/* FAST X
  for(i=0;i<*WRITE_STEP;i++){
    dispArray[i] = i;
    dispArray[i] *= sizeof(float)*6*sgt_numsta;
  }
  MPI_Type_create_hindexed(*WRITE_STEP, ones, dispArray, filetype, &filetype);
  MPI_Type_commit(&filetype);
*/
    free(ones);
    free(dispArray);

    // BEWARE that if there are other threads requesting other receivers
    // to be transferred to host, numsta is superset of sgt_numsta
    *numsta = sgt_numsta;
    *numsta_local = sgt_numsta_local;
    //*sta = &sgt_sta;
    //*sta_indices = &sgt_sta_indices;
    *areRegularReceivers = 0;
    sgt_buf = (float*)malloc(sizeof(float)*sgt_numsta_local*6*(*WRITE_STEP));

    printf("%d) In sgt, after pmcl3d's parameters are set\n",*rank);
    printf("%d) In sgt, sta:%d,%d,%d; %d,%d,%d; ...; %d,%d,%d\n",
           *rank, (*sta)[0],(*sta)[1],(*sta)[2],(*sta)[3],(*sta)[4],(*sta)[5],
           (*sta)[*numsta_local-3],(*sta)[*numsta_local-2],(*sta)[*numsta_local-1]);

    if(sgt_igreen != -1){
        sgt_c1 = (float*)malloc(sizeof(float)*sgt_numsta_local);
        sgt_c2 = (float*)malloc(sizeof(float)*sgt_numsta_local);
        sgt_mu = (float*)malloc(sizeof(float)*sgt_numsta_local);
        for(i=0;i<sgt_numsta_local;i++){
            ind0 = (*sta)[i*3]   +2+4*loop;
            ind1 = (*sta)[i*3+1] +2+4*loop;
            ind2 = (*sta)[i*3+2] +2+4*loop;
            tmplam = lam[ind0][ind1][ind2];
            tmpmu = mu[ind0][ind1][ind2];
            sgt_mu[i] = tmpmu;
            sgt_c1[i] = (1./tmplam+1./tmpmu)/(1./tmpmu*(3./tmplam+2./tmpmu));
            sgt_c2[i] = -1./tmplam/(2./tmpmu*(3./tmplam+2./tmpmu));
        }
    }

    (*sgt_args)->rank = rank;
    (*sgt_args)->cur_step_ptr = cur_step_ptr;
    (*sgt_args)->nt = nt;
    (*sgt_args)->NTISKP = NTISKP;
    (*sgt_args)->WRITE_STEP = WRITE_STEP;
    (*sgt_args)->ntiskp_mutex = ntiskp_mutex;
    (*sgt_args)->veldata_rwlock = veldata_rwlock;
    (*sgt_args)->ntiskp_cond = ntiskp_cond;
    (*sgt_args)->stream_IO = stream_IO;
    (*sgt_args)->tmpbuf = tmpbuf;
    (*sgt_args)->fh = fh;
    (*sgt_args)->displacement = displacement;
    (*sgt_args)->filetype = filetype;

    return;
}