#include <iostream>
//#include <MVMultiply.h>
#include <cmath>

using namespace std;

double* mvMultiply(double** matrix, double* vector, int size);
void printMatrix(double** matrix, int size);
void printVector(double* vector, int size);

int main()
{
    int size = 5;
    double** matrix = new double*[size];
    double** diagmatrix = new double*[size];
    double* vector = new double[size];

    for (int i = 0; i < size; i++)
    {
        matrix[i] = new double[size];
        diagmatrix[i] = new double[size];
    }

    printMatrix(diagmatrix, size);
    printVector(vector, size);
    

    double* solVect = mvMultiply(matrix, vector, size);
    printVector(solVect, size);
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