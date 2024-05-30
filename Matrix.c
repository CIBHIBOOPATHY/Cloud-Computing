/* This program performs simple matrix operations in a parallel distributed environment */
#include<stdio.h>
#include<mpi.h>

//Define number of rows of the square matrices
#define N 1000

/* Define matrices
 *
 * A,B - Input matrices
 * C - Matrix addition result matrix
 * D - Matrix multiplication result matrix
 * valC - Validation matrix for matrix addition
 * valD - Validation matrix for matrix multiplication
 */

double A[N][N], B[N][N], C[N][N], D[N][N], valC[N][N], valD[N][N];

int i,j,k;

/*
 * Function: mat_gen
 * -----------------
 *   generates matrices for matrix operations
 *
 *   returns: None
 */

void mat_gen(){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            A[i][j]=1+j;
            B[i][j]=1;
        }
    }
}

/*
 * Function: matrixop
 * ------------------
 *   Performs matrix operations on matrices A and B - Addition and Multiplication.
 *
 *   lower_index: Starting row index of matrices
 *   upper_index: (Index of last row + 1)
 *   val: validation number to perform validation matrix operation or not
 *
 *   returns: None
 */
void matrixop(int lower_index, int upper_index, int val){
    for (i = lower_index; i < upper_index; i++){
        for (j = 0; j < N; j++) {
            if(val==1){
                valC[i][j]=A[i][j]+B[i][j];
            }
            else{
                C[i][j]=A[i][j]+B[i][j];
            }
            for (k = 0; k < N; k++) {
                if(val==1){
                    valD[i][j] += A[i][k] * B[k][j];
                }
                else{
                    D[i][j] += A[i][k] * B[k][j];
                }
            }
        }
    }
}

/*
 * Function: matrixval
 * -------------------
 *   Performs validation of results obtained from matrix operations by comparing each element
 *   of validation matrix and result matrix
 *
 *   returns: None
 */

void matrixval(){
    int val_add=0;
    int val_mul=0;

    printf("\nValidating for matrix operations with a single node matrix operation results...\n");
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            if(valC[i][j]==C[i][j]){
                val_add+=1;
            }
            if(valD[i][j]==D[i][j]){
                val_mul+=1;
            }
        }
    }
    if(val_add==(N*N)){
        printf("\nMatrix addition is validated...\n");
    }
    if(val_mul==(N*N)){
        printf("\nMatrix multiplication is also validated...\n");
    }
    if(val_mul!=(N*N) || val_add!=(N*N)){
        printf("\nValidation failed...\n");
    }
}

/*
 * Function: printmatrix
 * ---------------------
 *   Prints input matrices, result matrices and validation matrices
 *
 *   returns: None
 */

void printmatrix(){
    for(int l=1; l<7;l++){
        if(l==1){
            printf("\nMatrix A\n");
        }
        if(l==2){
            printf("\nMatrix B\n");
        }
        if(l==3){
            printf("\nMatrix Addition Result A+B\n");
        }
        if(l==4){
            printf("\nValidation Matrix Addition(single node) result A+B\n");
        }
        if(l==5){
            printf("\nMatrix Multiplication result A*B\n");
        }
        if(l==6){
            printf("\nValidation Matrix Multiplication(single node) result A*B\n");
        }
        for ( i = 0; i < N; i++){
            printf("\n");
            for (int j = 0; j < N; j++){
                if(l==1){
                    printf("%.2f  ", A[i][j]);
                }
                if(l==2){
                    printf("%.2f  ", B[i][j]);
                }
                if(l==3){
                    printf("%.2f  ", C[i][j]);
                }
                if(l==4){
                    printf("%.2f  ", valC[i][j]);
                }
                if(l==5){
                    printf("%.2f  ", D[i][j]);
                }
                if(l==6){
                    printf("%.2f  ", valD[i][j]);
                }
            }
        }
        printf("\n\n");
    }
}

int main(int argc, char *argv[])
{
    int rank,np; // defines rank and number of process
    double begin_time, end_time; // start time and end time for matrix operations

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Status status;
    MPI_Request request;

    // matrix generation
    mat_gen();

    // validation matrix calculation
    matrixop(0,N,1);

    begin_time=MPI_Wtime();

    // broadcast matrix B to all nodes
    MPI_Bcast(&B, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // define the work package size for matrix decomposition
    int work_size=N/np;
    int lower_index=0;
    int upper_index=lower_index+work_size;
    /*---------------------------------------------------MASTER-NODE---------------------------------------------------------*/
    if(rank==0){
        // matrix operations on first work package
        matrixop(lower_index,upper_index,0);

        printf("\nMatrix is decomposed and sent to multiple processing nodes\n");

        // sending remaining work packages to different nodes
        for(i=1;i<np;i++){
            lower_index=i*work_size;
            if((i+1)==np){
                upper_index=N;
            }
            else{
                upper_index=lower_index+work_size;
            }

            MPI_Isend(&lower_index, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
            MPI_Isend(&upper_index, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &request);
            MPI_Isend(&A[lower_index][0], (upper_index - lower_index) * N, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &request);
        }
    }
	/*---------------------------------------------------WORKER-NODE---------------------------------------------------------*/
    if(rank>0){
        // receive the decomposed matrix work packages
        MPI_Recv(&lower_index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&upper_index, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&A[lower_index][0], (upper_index - lower_index) * N, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);

        matrixop(lower_index,upper_index,0); // perform matrix operation on work packages

        // send back the work packages after performing matrix operations to master node
        MPI_Isend(&lower_index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
        MPI_Isend(&upper_index, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &request);
        MPI_Isend(&C[lower_index][0], (upper_index - lower_index) * N, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &request);
        MPI_Isend(&D[lower_index][0], (upper_index - lower_index) * N, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &request);
    }

    /*---------------------------------------------------MASTER-NODE---------------------------------------------------------*/
    if (rank==0){
        printf("\nGathering matrix operation results from worker nodes...\n");
        for (i = 1; i < np; i++){
            // gather the result work packages
            MPI_Recv(&lower_index, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&upper_index, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&C[lower_index][0], (upper_index - lower_index) * N, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&D[lower_index][0], (upper_index - lower_index) * N, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);
        }
        end_time = MPI_Wtime();
        printf("\nRun Time for performing simple matrix operations = %f seconds\n", end_time - begin_time); // calculate t$

        printmatrix(); // print the matrices
        matrixval(); // check validation for results obtained
    }
    MPI_Finalize(); // ending the MPI operations
    return 0;
}
