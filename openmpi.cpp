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


/*
 * Procedure merge_array
 *
 * Merges two subarrays of n with n[left..mid] and n[mid+1..right]
 * to n itself, with n now ordered ascendingly
 * */
void merge_array(int *n, int left, int mid, int right) {
	int n_left = mid - left + 1;
	int n_right = right - mid;
	int iter_left = 0, iter_right = 0, iter_merged = left;
	int arr_left[n_left], arr_right[n_right];

	for (int i = 0; i < n_left; i++) {
		arr_left[i] = n[i + left];
	}

	for (int i = 0; i < n_right; i++) {
		arr_right[i] = n[i + mid + 1];
	}

	while (iter_left < n_left && iter_right < n_right) {
		if (arr_left[iter_left] <= arr_right[iter_right]) {
			n[iter_merged] = arr_left[iter_left++];
		} else {
			n[iter_merged] = arr_right[iter_right++];
		}
		iter_merged++;
	}

	while (iter_left < n_left)  {
		n[iter_merged++] = arr_left[iter_left++];
	}
	while (iter_right < n_right) {
		n[iter_merged++] = arr_right[iter_right++];
	} 
}

/* 
 * Procedure merge_sort
 *
 * Sorts array n with merge sort algorithm
 * */
void merge_sort(int *n, int left, int right) {
	if (left < right) {
		int mid = left + (right - left) / 2;

		merge_sort(n, left, mid);
		merge_sort(n, mid + 1, right);

		merge_array(n, left, mid, right);
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
    int* myArray;
    int sizeMyArray=0;

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

        //TODO initate the value of myArray and sizeMyArray with the result of convolution. Array with one element i think.
    }

    /* Sorting Logic Starts From Here!1!1!1! */
    
    //variable to help
    int divisor = 2;
    int divisor_difference = 1;

    //start process
    while (divisor <= size){
        if (rank % divisor == 0){
            int rankPartner = rank + divisor_difference
            if (rankPartner < size){    //check condition is the partner exist
                //receive size of array
                int sizeOfRecvArray;
                MPI_Recv(&sizeOfRecvArray, 1, MPI_INT, rankPartner, 1, MPI_COMM_WORLD, 0);  //tag 1 only for differentiate this message to get size of array
                
                //receive array
                int arrayOfInt[sizeOfRecvArray];
                MPI_Recv(arrayOfInt, sizeOfRecvArray, MPI_INT, rankPartner, 2, MPI_COMM_WORLD, 0); //tag 2 only for differentiate this message to get the array

                //unite myArray with received array
                sizeMyArray = sizeMyArray + sizeOfRecvArray
                int* tempArray;
                tempArray = (int *)malloc((sizeMyArray)*sizeof(int));
                for (int i=0; i<sizeMyArray; i++){
                    tempArray[i] = myArray[i];
                }
                for (int i=0; i<sizeOfRecvArray; i++){
                    tempArray[i+sizeMyArray] = arrayOfInt[i];
                }
                free(myArray);
                myArray = tempArray;
                tempArray = NULL;

                //sort array
                merge_sort(myArray, 0, sizeMyArray-1);
            }
        }else{
            int rankPartner = rank - divisor_difference;

            //send size of array
            MPI_Send(&sizeMyArray, 1, MPI_INT, rankPartner, 1, MPI_COMM_WORLD);

            //send the array
            MPI_Send(myArray, sizeMyArray, MPI_INT, rankPartner, 2, MPI_COMM_WORLD);
        }

        divisor *= 2
        divisor_difference *= 2
    }

    MPI_Finalize();
}