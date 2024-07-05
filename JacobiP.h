#include <iostream>
#include <iomanip>
#include <SCmathlib.h>
#include <mpi.h>

using namespace std;

const int matsize = 9; // N = 9: h = 0.1

double* mvMultiplyP(double** matrix, double* vector, int row, int col, int rank, int nodesize);
void JacobiP(double** matrix, double* vector, double* solvect, int matsize, int rank, int nodesize, double TOL, int maxIter);
void printMatrix(double** matrix, int size);
void printVector(double* vector, int size);

void JacobiP(double** matrix, double* vector, double* solvect, int matsize, int rank, int nodesize, double TOL, int maxIter)
{
    /*             Jacobi Algorithm:              */
    /* Solve for x such that:                     */
    /* Dx(i+1) = b - (A-D)x(i)                    */
    /* where x(i) is the initial guess,           */
    /*       D is the diagonal matrix of Matrix A */
    /*       b is the vector such that Ax = b     */
    /*       x(i+1) is the solution vector        */

    double** diagmat = new double*[matsize];
    double** invdiag = new double*[matsize];

    double* truesol = new double[matsize];

    for (int i = 0; i < matsize; i++)
    {
        diagmat[i] = new double[matsize];
        invdiag[i] = new double[matsize];

        for (int j = 0; j < matsize; j++)
        {
            if (i == j)
            {
                diagmat[i][j] = matrix[i][j];
                invdiag[i][j] = 1.0/diagmat[i][j];
            }

            else
            {
                diagmat[i][j] = 0;
                invdiag[i][j] = 0;
            }
        }
    }

    // Jacobi iteration
    int iter;

    for (iter = 0; iter < maxIter; iter++)
    {
        double** dif1 = new double*[matsize];

        for (int i = 0; i < matsize; i++)
        {
            dif1[i] = new double[matsize];
            for (int j = 0; j < matsize; j++)
            {
                dif1[i][j] = -diagmat[i][j] + matrix[i][j];
            }
        }

        double* prod1 = mvMultiplyP(dif1, solvect, matsize, matsize, rank, nodesize);

        double* dif2 = new double[matsize];

        for (int i = 0; i < matsize; i++)
        {
            dif2[i] = vector[i] - prod1[i];
        }

        truesol = mvMultiplyP(invdiag, dif2, matsize, matsize, rank, nodesize);

        double sumerr = 0.0;

        for (int i = 0; i < matsize; i++)
        {
            sumerr += pow(truesol[i] - solvect[i],2);
        }

        if (sqrt(sumerr) < TOL)
            break;

        // Prepare for next iteration
        for (int i = 0; i < matsize; i++)
        {
            solvect[i] = truesol[i];
            delete[] dif1[i];
        }

        delete[] dif1;
        delete[] dif2;
        delete[] prod1;
    }

    if (rank == 0)
    {
        cout << "Solution ";
        printVector(solvect, matsize);
        cout << "Achieved in: " << iter << " iterations." << endl;
    }

    for (int i = 0; i < matsize; i++)
    {
        delete[] diagmat[i];
        delete[] invdiag[i];
    }

    delete[] diagmat;
    delete[] invdiag;
    delete[] truesol;
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

    MPI_Bcast(vector, col, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Local matrix-vector multiplication
    for (int i = 0; i < localRowSize; i++)
    {
        prodVectLocal[i] = 0;

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
            if (i == nodesize - 1)
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
