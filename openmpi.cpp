#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define NMAX 100
#define MATRIXNSEND 10002
#define DATAMAX 1000
#define DATAMIN -1000

/*
 * Struct Matrix
 *
 * Matrix representation consists of matrix data
 * and effective dimensions
 * */
typedef struct Matrix
{
    int mat[NMAX][NMAX]; // Matrix cells
    int row_eff;         // Matrix effective row
    int col_eff;         // Matrix effective column
} Matrix;

/*
 * Procedure init_matrix
 *
 * Initializing newly allocated matrix
 * Setting all data to 0 and effective dimensions according
 * to nrow and ncol
 * */
void init_matrix(Matrix *m, int nrow, int ncol)
{
    m->row_eff = nrow;
    m->col_eff = ncol;

    for (int i = 0; i < m->row_eff; i++)
    {
        for (int j = 0; j < m->col_eff; j++)
        {
            m->mat[i][j] = 0;
        }
    }
}

/*
 * Function input_matrix
 *
 * Returns a matrix with values from stdin input
 * */
Matrix input_matrix(int nrow, int ncol)
{
    Matrix input;
    init_matrix(&input, nrow, ncol);

    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < ncol; j++)
        {
            scanf("%d", &input.mat[i][j]);
        }
    }

    return input;
}

/*
 * Procedure print_matrix
 *
 * Print matrix data
 * */
void print_matrix(Matrix *m)
{
    for (int i = 0; i < m->row_eff; i++)
    {
        for (int j = 0; j < m->col_eff; j++)
        {
            printf("%d ", m->mat[i][j]);
        }
        printf("\n");
    }
}

int temp[MATRIXNSEND];
int pivot = 1;

void compose_matrix(Matrix *m)
{
    m->row_eff = temp[0];
    m->col_eff = temp[1];
    for (int i = 0; i < row_eff; i++)
    {
        for (int j = 0; j < col_eff; j++)
        {
            m->mat[i][j] = temp[pivot + (i * m->row_eff) + j + 1];
        }
    }
}

void decompose_matrix(Matrix *m)
{
    for (int i = 0; i < MATRIXNSEND; i++)
    {
        temp[i] = 0;
    }
    temp[0] = m->row_eff;
    temp[1] = m->col_eff;

    for (int i = 0; i < row_eff; i++)
    {
        for (int j = 0; j < col_eff; j++)
        {
            temp[pivot + (i * m->row_eff) + j + 1] = m->mat[i][j];
        }
    }
}

int main(int argc, char *argv[])
{
    int size;
    int rank;
    Matrix inputMatrix;
    Matrix *kernelMatrix;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0)
    {
        int inputRow, inputCol;
        scanf("%d %d", &inputRow, &inputCol);
        inputMatrix = input_matrix(inputRow, inputCol);

        int numTargets, targetRow, targetCol;
        scanf("%d %d %d", &numTargets, &targetRow, &targetCol);
        kernelMatrix = (Matrix *)malloc(numTargets * sizeof(Matrix));

        for (int i = 0; i < numTargets; i++)
        {
            kernelMatrix[i] = input_matrix(targetRow, targetCol);
        }

        // Bagian bagi-bagi kerjaan
        int targetMatrixPerRank = numTargets / size;
        int forEachSlave = targetMatrixPerRank;
        int forMaster = numTargets - (targetMatrixPerRank * (size - 1));
        // Jadi numTargets bakal dibagi bersama size
        // Trus kerjaan per slave itu floornya, sedangkan kerjaan master sisanya

        for (int i = 1; i < size; i++)
        {
            decompose_matrix(&inputMatrix);
            // kasih input matrix ke slavenya
            MPI_Send(temp, MATRIXNSEND, MPI_INT, i, 0, MPI_COMM_WORLD);
            // kasih tau jumlah kerjaan ke slavenya
            MPI_Send(&forEachSlave, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

            for (int j = 0; j < forEachSlave; j++)
            {
                decompose_matrix(&(kernelMatrix[j + (i - 1) * forEachSlave + forMaster]));
                // kasih tau matrix kernel kerjaan ke j ke slavenya
                MPI_Send(temp, MATRIXNSEND, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }

        for (int i = 0; i < forMaster; i++)
        {
            // Do kerjaan master
            // KernelMatrix buat master itu kernelMatrix[0] sampai kernelMatrix[forMaster-1] atau sejumlah forMaster
        }
    }
    else
    {
        // Dapetin matrix input dari master
        MPI_Recv(temp, MATRIXNSEND, MPI_INT, 0, 0, MPI_COMM_WORLD, 0);
        compose_matrix(&inputMatrix);

        int numTargets;
        // minta jumlah kerjaan dari master
        MPI_Recv(&numTargets, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, 0);
        kernelMatrix = (Matrix *)malloc(numTargets * sizeof(Matrix));

        for (int i = 0; i < numTargets; i++)
        {
            // minta matrix kernel dari master
            MPI_Recv(temp, MATRIXNSEND, MPI_INT, 0, 0, MPI_COMM_WORLD, 0);
            compose_matrix(&(kernelMatrix[i]));
        }

        for (int i = 0; i < numTargets; i++)
        {
            // Do kerjaan slave
            // KernelMatrix buat slave ini itu kernelMatrix[0] sampai kernelMatrix[numTargets-1] atau sejumlah numTargets
        }
    }

    MPI_Finalize();
}