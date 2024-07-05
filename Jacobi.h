#include <iostream>
#include <cmath>
#include <chrono>

using namespace std;
using namespace std::chrono;

double* mvMultiply(double** matrix, double* vector, int row, int col);
void Jacobi(double** matrix, double* vector, double* solvect, int matsize, double TOL, int maxIter);
void printMatrix(double** matrix, int size);
void printVector(double* vector, int size);

void Jacobi(double** matrix, double* vector, double* solvect, int matsize, double TOL, int maxIter)
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

    double* truesol;

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

        double* prod1 = mvMultiply(dif1, solvect, matsize, matsize);

        double* dif2 = new double[matsize];

        for (int i = 0; i < matsize; i++)
        {
            dif2[i] = vector[i] - prod1[i];
        }

        truesol = mvMultiply(invdiag, dif2, matsize, matsize);

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

    cout << "Solution ";
    printVector(solvect, matsize);
    cout << "Achieved in: " << iter << " iterations." << endl;

    for (int i = 0; i < matsize; i++)
    {
        delete[] diagmat[i];
        delete[] invdiag[i];
    }

    delete[] diagmat;
    delete[] invdiag;
    delete[] truesol;
}

double* mvMultiply(double** matrix, double* vector, int row, int col)
{
    double* solVect = new double[row];

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
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
