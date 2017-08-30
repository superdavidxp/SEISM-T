/*
 *
 */

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "seism_t.h"
#include "math.h"
#include "pthread.h"

#define BLOCK_SIZE_X 2
#define BLOCK_SIZE_Y 2
#define BLOCK_SIZE_Z 128

void calcRegularReceiverPoints(int *rec_nbgx, int *rec_nedx, int *rec_nbgy, int *rec_nedy, int *rec_nbgz, int *rec_nedz,
                               int *rec_nxt, int *rec_nyt, int *rec_nzt, MPI_Offset *displacement,
                               long int nxt, long int nyt, long int nzt, int rec_NX, int rec_NY, int rec_NZ,
                               int NBGX, int NEDX, int NSKPX, int NBGY, int NEDY, int NSKPY,
                               int NBGZ, int NEDZ, int NSKPZ, int *coord);

void calculate_IO_BLOCK_SIZE(int *IO_BLOCK_SIZE_X, int *IO_BLOCK_SIZE_Y, int *IO_BLOCK_SIZE_Z,
                             int rec_nxt, int rec_nyt, int rec_nzt);

Grid3D Alloc3D(int nx, int ny, int nz);

Grid1D Alloc1D(int nx);

int main(int argc, char** argv)
{
    int i        = -1;
    int j        = -1;
    int k        = -1;
    int tmpSize  = -1;
    int err      =  0;
    int irank    = -1;
    int nrank    = -1;
    int provided = -1;
    int cur_step = -1;
    int idtmp    = -1;
    int tmpInd   = -1;
    float *tmpbuf;

    MPI_Comm MCW, MC1;
    MPI_Offset displacement;
    MPI_Datatype filetype;
    MPI_File fh;

    seism_t_mpi_init(argc, argv, &provided);
    err = MPI_Comm_size(MPI_COMM_WORLD, &nrank);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &irank);

    err = MPI_Comm_dup(MPI_COMM_WORLD,&MCW);
    err = MPI_Barrier(MCW);

    printf("RANK : %d, of %d : provided = %d\n", irank, nrank, provided);

    pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
    pthread_mutex_init(&mpi_mutex,NULL);
    pthread_cond_init(&mpi_cond,NULL);
    pthread_mutex_init(&ntiskp_mutex,NULL);
    pthread_cond_init(&ntiskp_cond,NULL);
    pthread_mutex_init(&vel_out_mutex,NULL);
    pthread_cond_init(&vel_out_cond,NULL);
    pthread_rwlock_init(&veldata_rwlock,NULL);

    float TMAX       = 0.4;
    float DH         = 40.0;
    float DT         = 0.02;
    float ARBC       = 0.95;
    float PHT        = 0.1;
    int   NPC        = 0;
    int   ND         = 10;
    int   NSRC       = 1;
    int   NST        = 101;
    int   NVAR       = 3;
    int   NVE        = 1;
    int   MEDIASTART = 0;
    int   IFAULT     = 1;
    int   READ_STEP  = 10;
    int   NTISKP     = 1;
    int   WRITE_STEP = 5;
    int   NX         = 128;
    int   NY         = 128;
    int   NZ         = 128;
    int   PX         = 1;
    int   PY         = 2;
    int   NBGX       = 1;
    int   NEDX       = 128;
    int   NSKPX      = 1;
    int   NBGY       = 1;
    int   NEDY       = 128;
    int   NSKPY      = 1;
    int   NBGZ       = 1;
    int   NEDZ       = 128;
    int   NSKPZ      = 1;
    float FAC        = 1.0;
    float Q0         = 150.0;
    float EX         = 0.0;
    float FP         = 1.0;
    int   IDYNA      = 0;
    int   SoCalQ     = 1;
    char  INSRC[512];
    char  INVEL[512];
    char  CHKFILE[512];

    int dim[2], period[2], coord[2], reorder;

    int nxt   = NX/PX;
    int nyt   = NY/PY;
    int nzt   = NZ;
    int nt    = (int)(TMAX/DT) + 1;
    dim[0]    = PX;
    dim[1]    = PY;
    period[0] = 0;
    period[1] = 0;
    reorder   = 1;

    int x_rank_L = -1;
    int x_rank_R = -1;
    int y_rank_F = -1;
    int y_rank_B = -1;

    err = MPI_Cart_create(MCW, 2, dim, period, reorder, &MC1);
    err = MPI_Cart_shift(MC1, 0,  1,  &x_rank_L, &x_rank_R );
    err = MPI_Cart_shift(MC1, 1,  1,  &y_rank_F, &y_rank_B );
    err = MPI_Cart_coords(MC1, irank, 2, coord);
    err = MPI_Barrier(MCW);

    int xls, xre, xvs, xve, xss1, xse1, xss2, xse2, xss3, xse3;
    int yfs, yfe, ybs, ybe, yls, yre;

    if(x_rank_L<0)
        xls = 2+4*loop;
    else
        xls = 4*loop;

    if(x_rank_R<0)
        xre = nxt+4*loop+1;
    else
        xre = nxt+4*loop+3;

    xvs   = 2+4*loop;
    xve   = nxt+4*loop+1;

    xss1  = xls;
    xse1  = 4*loop+3;
    xss2  = 4*loop+4;
    xse2  = nxt+4*loop-1;
    xss3  = nxt+4*loop;
    xse3  = xre;

    if(y_rank_F<0)
        yls = 2+4*loop;
    else
        yls = 4*loop;

    if(y_rank_B<0)
        yre = nyt+4*loop+1;
    else
        yre = nyt+4*loop+3;

    yfs  = 2+4*loop;
    yfe  = 2+8*loop-1;
    ybs  = nyt+2;
    ybe  = nyt+4*loop+1;

    // same for each processor:
    if (NEDX==-1) NEDX = NX;
    if (NEDY==-1) NEDY = NY;
    if (NEDZ==-1) NEDZ = NZ;
    NEDX = NEDX-(NEDX-NBGX)%NSKPX;
    NEDY = NEDY-(NEDY-NBGY)%NSKPY;
    NEDZ = NEDZ-(NEDZ-NBGZ)%NSKPZ;

    int rec_NX;
    int rec_NY;
    int rec_NZ;
    int rec_nxt;
    int rec_nyt;
    int rec_nzt;
    int rec_nbgx;   // 0-based indexing, however NBG* is 1-based
    int rec_nedx;   // 0-based indexing, however NED* is 1-based
    int rec_nbgy;   // 0-based indexing
    int rec_nedy;   // 0-based indexing
    int rec_nbgz;   // 0-based indexing
    int rec_nedz;   // 0-based indexing
    int IO_BLOCK_SIZE_X;
    int IO_BLOCK_SIZE_Y;
    int IO_BLOCK_SIZE_Z;

    // number of receiving points in total
    rec_NX = (NEDX-NBGX)/NSKPX+1;
    rec_NY = (NEDY-NBGY)/NSKPY+1;
    rec_NZ = (NEDZ-NBGZ)/NSKPZ+1;

    calcRegularReceiverPoints(&rec_nbgx, &rec_nedx,
                              &rec_nbgy, &rec_nedy, &rec_nbgz, &rec_nedz,
                              &rec_nxt, &rec_nyt, &rec_nzt, &displacement,
                              (long int)nxt,(long int)nyt,(long int)nzt,
                              rec_NX, rec_NY, rec_NZ,
                              NBGX,NEDX,NSKPX, NBGY,NEDY,NSKPY, NBGZ,NEDZ,NSKPZ, coord);

    printf("%d = (%d,%d)) NX,NY,NZ=%d,%d,%d\nnxt,nyt,nzt=%d,%d,%d\nrec_N=(%d,%d,%d)\nrec_nxt,=%d,%d,%d\nNBGX,SKP,END=(%d:%d:%d),(%d:%d:%d),(%d:%d:%d)\nrec_nbg,ed=(%d,%d),(%d,%d),(%d,%d)\ndisp=%ld\n",
           irank,coord[0],coord[1],NX,NY,NZ,nxt,nyt,nzt,
           rec_NX, rec_NY, rec_NZ, rec_nxt, rec_nyt, rec_nzt,
           NBGX,NSKPX,NEDX,NBGY,NSKPY,NEDY,NBGZ,NSKPZ,NEDZ,
           rec_nbgx,rec_nedx,rec_nbgy,rec_nedy,rec_nbgz,rec_nedz,(long int)displacement);

    // calculate IO_BLOCK_SIZE_X,Y,Z
    calculate_IO_BLOCK_SIZE(&IO_BLOCK_SIZE_X, &IO_BLOCK_SIZE_Y, &IO_BLOCK_SIZE_Z,
                            rec_nxt, rec_nyt, rec_nzt);

    int maxNX_NY_NZ_WS = -1;
    maxNX_NY_NZ_WS = (rec_NX>rec_NY?rec_NX:rec_NY);
    maxNX_NY_NZ_WS = (maxNX_NY_NZ_WS>rec_NZ?maxNX_NY_NZ_WS:rec_NZ);
    maxNX_NY_NZ_WS = (maxNX_NY_NZ_WS>WRITE_STEP?maxNX_NY_NZ_WS:WRITE_STEP);

    printf("|    maxNX_NY_NZ_WS = %d\n", maxNX_NY_NZ_WS);

    int ones[maxNX_NY_NZ_WS];
    MPI_Aint dispArray[maxNX_NY_NZ_WS];
    for (i=0;i<maxNX_NY_NZ_WS;++i) {
        ones[i] = 1;
    }

    err = MPI_Type_contiguous(rec_nxt, MPI_FLOAT, &filetype);
    err = MPI_Type_commit(&filetype);
    for (i=0;i<rec_nyt;i++) {
        dispArray[i] = sizeof(float);
        dispArray[i] = dispArray[i]*rec_NX*i;
    }
    err = MPI_Type_create_hindexed(rec_nyt, ones, dispArray, filetype, &filetype);
    err = MPI_Type_commit(&filetype);
    for(i=0;i<rec_nzt;i++){
        dispArray[i] = sizeof(float);
        dispArray[i] = dispArray[i]*rec_NY*rec_NX*i;
    }
    err = MPI_Type_create_hindexed(rec_nzt, ones, dispArray, filetype, &filetype);
    err = MPI_Type_commit(&filetype);
    for(i=0;i<WRITE_STEP;i++){
        dispArray[i] = sizeof(float);
        dispArray[i] = dispArray[i]*rec_NZ*rec_NY*rec_NX*i;
    }
    err = MPI_Type_create_hindexed(WRITE_STEP, ones, dispArray, filetype, &filetype);
    err = MPI_Type_commit(&filetype);
    MPI_Type_size(filetype, &tmpSize);
    printf("filetype size (supposedly = rec_nxt*nyt*nzt*WS*4 = %d) = %d\n", (int)(rec_nxt*rec_nyt*rec_nzt*WRITE_STEP*4), tmpSize);
    printf("Allocate device IO buffers.\n");

    Grid1D Bufx=NULL, Bufy=NULL, Bufz=NULL;

    Grid3D u1=NULL, v1=NULL, w1=NULL;
    u1  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*align);
    v1  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*align);
    w1  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*align);

    float *d_u1, *d_v1, *d_w1;
    if (irank==0) printf("Allocate device velocity and stress pointers and copy.\n");
    long int num_bytes = sizeof(float)*(nxt+4+8*loop)*(nyt+4+8*loop)*(nzt+2*align);

    num_bytes = rec_nxt*rec_nyt*rec_nzt;
    Bufx  = Alloc1D(num_bytes*WRITE_STEP);
    Bufy  = Alloc1D(num_bytes*WRITE_STEP);
    Bufz  = Alloc1D(num_bytes*WRITE_STEP);
    num_bytes = rec_nxt*rec_nyt*rec_nzt*3;

    tmpbuf = (float*) malloc (sizeof(float)*num_bytes);

    printf("RANK : %d, IO_BLOCK_SIZE = %d,%d,%d\n", irank, IO_BLOCK_SIZE_X, IO_BLOCK_SIZE_Y, IO_BLOCK_SIZE_Z);

    int cur_step_ntiskp = 0;
    char filename[50];
    char filenamebasex[25] = "./SX";
    char filenamebasey[25] = "./SY";
    char filenamebasez[25] = "./SZ";

    //MPI_Comm_split(MPI_COMM_WORLD, velout_MPIcolor, *rank, &velout_MCwriters);
    MPI_Comm_dup(MPI_COMM_WORLD, &velout_MCwriters);
    
    vel_out_init(&velout_args, &irank, Bufx, Bufy, Bufz,
                 &rec_nxt, &rec_nyt, &rec_nzt,
                 &nt, &NTISKP, &WRITE_STEP,
                 &cur_step, &fh,
                 &displacement, &filetype,
                 filename, filenamebasex,
                 filenamebasey, filenamebasez);

    int ioutput = 0;

    //  Main Loop Starts
    if (NPC==0 && NVE==1)
    {
        for (cur_step=1;cur_step<=nt;cur_step++)
        {
            if (irank == 0) printf("|    Time Step =                   %7d  OF  Total Timesteps = %d\n", cur_step, nt);

            // calculation
            sleep(0.05);

            // this whole loop go in to pthread library
            if ( fmod(cur_step, NTISKP) == 0 )
            {
                num_bytes = sizeof(float)*(nxt+4+8*loop)*(nyt+4+8*loop)*(nzt+2*align);
                idtmp = ((cur_step/NTISKP+WRITE_STEP-1)%WRITE_STEP);
                idtmp = idtmp*rec_nxt*rec_nyt*rec_nzt;
                tmpInd = idtmp;

                for(k=nzt+align-1 - rec_nbgz; k>=nzt+align-1 - rec_nedz; k=k-NSKPZ)
                {
                    for(j=2+4*loop + rec_nbgy; j<=2+4*loop + rec_nedy; j=j+NSKPY)
                    {
                        for(i=2+4*loop + rec_nbgx; i<=2+4*loop + rec_nedx; i=i+NSKPX)
                        {
                            Bufx[tmpInd] = u1[i][j][k];
                            Bufy[tmpInd] = v1[i][j][k];
                            Bufz[tmpInd] = w1[i][j][k];
                            tmpInd++;
                        }
                    }
                }

            } // if ( fmod(cur_step, NTISKP) == 0 )

            vel_out_pthread((void*) velout_args);

        } // for ( cur_step=1;cur_step<=nt;cur_step++ )



    } // if ( NPC==0 && NVE==1 )

    MPI_Barrier(MC1);
    MPI_Finalize();
    return 0;
}

void calcRegularReceiverPoints(int *rec_nbgx, int *rec_nedx, int *rec_nbgy, int *rec_nedy, int *rec_nbgz, int *rec_nedz,
                               int *rec_nxt, int *rec_nyt, int *rec_nzt, MPI_Offset *displacement,
                               long int nxt, long int nyt, long int nzt, int rec_NX, int rec_NY, int rec_NZ,
                               int NBGX, int NEDX, int NSKPX, int NBGY, int NEDY, int NSKPY,
                               int NBGZ, int NEDZ, int NSKPZ, int *coord)
{
    *displacement = 0;

    if(NBGX > nxt*(coord[0]+1))     *rec_nxt = 0;
    else if(NEDX < nxt*coord[0]+1)  *rec_nxt = 0;
    else{
        if(nxt*coord[0] >= NBGX){
            *rec_nbgx = (nxt*coord[0]+NBGX-1)%NSKPX;
            *displacement += (nxt*coord[0]-NBGX)/NSKPX+1;
        }
        else
            *rec_nbgx = NBGX-nxt*coord[0]-1;  // since rec_nbgx is 0-based
        if(nxt*(coord[0]+1) <= NEDX)
            *rec_nedx = (nxt*(coord[0]+1)+NBGX-1)%NSKPX-NSKPX+nxt;
        else
            *rec_nedx = NEDX-nxt*coord[0]-1;
        *rec_nxt = (*rec_nedx-*rec_nbgx)/NSKPX+1;
    }

    if(NBGY > nyt*(coord[1]+1))     *rec_nyt = 0;
    else if(NEDY < nyt*coord[1]+1)  *rec_nyt = 0;
    else{
        if(nyt*coord[1] >= NBGY){
            *rec_nbgy = (nyt*coord[1]+NBGY-1)%NSKPY;
            *displacement += ((nyt*coord[1]-NBGY)/NSKPY+1)*rec_NX;
        }
        else
            *rec_nbgy = NBGY-nyt*coord[1]-1;  // since rec_nbgy is 0-based
        if(nyt*(coord[1]+1) <= NEDY)
            *rec_nedy = (nyt*(coord[1]+1)+NBGY-1)%NSKPY-NSKPY+nyt;
        else
            *rec_nedy = NEDY-nyt*coord[1]-1;
        *rec_nyt = (*rec_nedy-*rec_nbgy)/NSKPY+1;
    }

    if(NBGZ > nzt) *rec_nzt = 0;
    else{
        *rec_nbgz = NBGZ-1;  // since rec_nbgz is 0-based
        *rec_nedz = NEDZ-1;
        *rec_nzt = (*rec_nedz-*rec_nbgz)/NSKPZ+1;
    }

    if(*rec_nxt == 0 || *rec_nyt == 0 || *rec_nzt == 0){
        *rec_nxt = 0;
        *rec_nyt = 0;
        *rec_nzt = 0;
    }

    // displacement assumes NPZ=1!
    *displacement *= sizeof(float);

    return;
}

void calculate_IO_BLOCK_SIZE(int *IO_BLOCK_SIZE_X, int *IO_BLOCK_SIZE_Y, int *IO_BLOCK_SIZE_Z,
                             int rec_nxt, int rec_nyt, int rec_nzt)
{
    int total_rec = rec_nxt*rec_nyt*rec_nzt;
    int BSIZES = BLOCK_SIZE_Z*BLOCK_SIZE_Y*BLOCK_SIZE_X;
    // default values: if(total_rec <= BSIZES)
    *IO_BLOCK_SIZE_X = 1;
    *IO_BLOCK_SIZE_Y = 1;
    *IO_BLOCK_SIZE_Z = 1;
    if(total_rec > BSIZES){
        float aggregation = ((float)total_rec)/BSIZES;
        if(aggregation <= rec_nzt)
            *IO_BLOCK_SIZE_Z = (int)aggregation+1;
        else{
            *IO_BLOCK_SIZE_Z = rec_nzt;
            aggregation = aggregation/rec_nzt;
            if(aggregation <= rec_nyt)
                *IO_BLOCK_SIZE_Y = (int)aggregation+1;
            else{
                *IO_BLOCK_SIZE_Y = rec_nyt;
                *IO_BLOCK_SIZE_X = (int)(aggregation/rec_nyt)+1;
            }
        }
    }
}

Grid3D Alloc3D(int nx, int ny, int nz)
{
    int i, j, k;
    Grid3D U = (Grid3D)malloc(sizeof(float**)*nx + sizeof(float *)*nx*ny +sizeof(float)*nx*ny*nz);

    if (!U){
        printf("Cannot allocate 3D float array\n");
        exit(-1);
    }
    for(i=0;i<nx;i++){
        U[i] = ((float**) U) + nx + i*ny;
    }

    float *Ustart = (float *) (U[nx-1] + ny);
    for(i=0;i<nx;i++)
        for(j=0;j<ny;j++)
            U[i][j] = Ustart + i*ny*nz + j*nz;

    for(i=0;i<nx;i++)
        for(j=0;j<ny;j++)
            for(k=0;k<nz;k++)
                U[i][j][k] = 0.0f;

    return U;
}

Grid1D Alloc1D(int nx)
{
    int i;
    Grid1D U = (Grid1D)malloc(sizeof(float)*nx);

    if (!U) {
        printf("Cannot allocate 2D float array\n");
        exit(-1);
    }

    for(i=0;i<nx;i++)
        U[i] = 0.0f;

    return U;
}
