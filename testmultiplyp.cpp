#include <iostream>
#include <cmath>
#include <mpi.h>

using namespace std;

double* mvMultiplyP(double** matrix, double* vector, int row, int col, int rank, int totalnodes);
void printMatrix(double** matrix, int size);
void printVector(double* vector, int size);

int main(int argc, char** argv)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int tprocess, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &tprocess);
    MPI_Comm_rank(comm, &rank);

    int size = 13;
    double** matrix = new double*[size];
    double** diagmatrix = new double*[size];
    double* vector = new double[size];

    for (int i = 0; i < size; i++)
    {
        matrix[i] = new double[size];
        diagmatrix[i] = new double[size];
    }

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrix[i][j] = i+j*2 + 1;
            vector[i] = floor(i*j*2) + 1;

            if (i == j)
                diagmatrix[i][j] = matrix[i][j];

            else
                diagmatrix[i][j] = 0;
        }
    }

    if (rank == 0)
    {
        printMatrix(matrix, size);
        printVector(vector, size);
    }
    
    double* solVect = mvMultiplyP(matrix, vector, size, size, rank, tprocess);
    if (rank == 0) printVector(solVect, size);

    for (int i = 0; i < size; i++)
    {
        delete[] matrix[i];
        delete[] diagmatrix[i];
    }

    delete[] matrix;
    delete[] diagmatrix;
    delete[] vector;
    delete[] solVect;

    MPI_Finalize();
}

double* mvMultiply(double** matrix, double* vector, int n)
{
    double* solVect = new double[n];

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            solVect[i] += matrix[i][j] * vector[j];
        }
    }

    return solVect;
}

void printMatrix(double** matrix, int size)
{
    cout << "Matrix: " << endl << endl;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            cout << matrix[i][j] << " ";
        }

        cout << endl << endl;;
    }
}

void printVector(double* vector, int size)
{
    cout << "Vector: " << endl << endl;
    for (int i = 0; i < size; i++)
    {
        cout << vector[i] << endl;
    }
    cout << endl;
}

double* mvMultiplyP(double** matrix, double* vector, int row, int col, int rank, int nodesize)
{
    MPI_Status status;

    int localNRow = row/nodesize;

    double* flatMatrix = new double[row*col];
    double* prodVect = new double[row];

    // reshape the matrix to be in one vector
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            flatMatrix[i*col + j] = matrix[i][j];
        }
    }

    /* distribute rows according to processors */
    
    int remainder = row%nodesize;

    int localRowSize;

    if (rank == nodesize - 1 && remainder != 0)
        localRowSize = localNRow + remainder;
    else
        localRowSize = localNRow;

    double* localMat = new double[localRowSize*col];
    double* prodVectLocal = new double[localRowSize];

    if (remainder == 0)
    {
        MPI_Scatter(flatMatrix, localNRow*col, MPI_DOUBLE, localMat, localNRow*col, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Case when unevenly-split rows by processors, and the last processor takes all
    else
    {
        int* sendcount = new int[nodesize];
        int* shiftcount = new int[nodesize];

        for (int i = 0; i < nodesize; i++)
        {
            if (i == nodesize - 1)
                sendcount[i] = (localNRow + remainder)*col;
            
            else
                sendcount[i] = (localNRow)*col;

            shiftcount[i] = localNRow*col * i; // by processor
        }

        MPI_Scatterv(flatMatrix, sendcount, shiftcount, MPI_DOUBLE, localMat, localRowSize*col, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        delete[] sendcount;
        delete[] shiftcount;
    }

    // Local matrix-vector multiplication
    for (int i = 0; i < localRowSize; i++)
    {
        for (int j = 0; j < col; j++)
        {
            prodVectLocal[i] += localMat[i*col + j] * vector[j];
        }
    } 

    if (remainder == 0)
        MPI_Allgather(prodVectLocal, localRowSize, MPI_DOUBLE, prodVect, localRowSize, MPI_DOUBLE, MPI_COMM_WORLD);

    else
    {
        int* recvcount = new int[nodesize];
        int* shiftcount = new int[nodesize];

        for (int i = 0; i < nodesize; i++)
        {
            if (i == nodesize - 1 && remainder != 0)
                recvcount[i] = localNRow + remainder;
            
            else
                recvcount[i] = localNRow;

            shiftcount[i] = i * localNRow;
        }

        MPI_Allgatherv(prodVectLocal, localRowSize, MPI_DOUBLE, prodVect, recvcount, shiftcount, MPI_DOUBLE, MPI_COMM_WORLD);
    
        delete[] recvcount;
        delete[] shiftcount;
    }

    delete[] flatMatrix;
    delete[] localMat;
    delete[] prodVectLocal;

    return prodVect;
}