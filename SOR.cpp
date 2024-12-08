#include <iostream>
#include <cmath>
#include <chrono>

using namespace std;
using namespace std::chrono;

const int matsize = 9;

void SOR(double** matrix, double* vector, double* initS, int n, double omega, double TOL, int maxIter);
void printMatrix(double** matrix, int size);
void printVector(double* vector, int size);

int main()
{
    double** matrix = new double*[matsize];
    double* vector = new double[matsize];

    double* initS = new double[matsize];

    double h = (1.0 - 0) / (matsize + 1); // (b-a)/(N+1)

    for (int i = 0; i < matsize; i++) 
    {
        matrix[i] = new double[matsize];

        for (int j = 0; j < matsize; j++) 
        {
            if (i == j) // diagonal part
                matrix[i][j] = 2.0;
            else if (j == i - 1 || j == i + 1) // shifted diagonal
                matrix[i][j] = -1.0;
            else // elsewhere
                matrix[i][j] = 0.0;
        }
    }

    double* xi = new double[matsize + 2];

    for (int i = 0; i < matsize + 2; i++) 
    {
        xi[i] = i * h;
    }

    for (int i = 0; i < matsize; i++) 
    {
        vector[i] = -(h * h) * sin(2.0 * M_PI * xi[i + 1]);
        initS[i] = 1.0; // initial guess
    }
    
    auto start = system_clock::now();
    SOR(matrix, vector, initS, matsize, 1.5, 1e-6, 10000);
    auto end = system_clock::now();

    duration<double> elapsedTime = end - start;
    cout << "Time elapsed: " << elapsedTime.count() << " seconds." << endl;
    
    for (int i = 0; i < matsize; i++)
        delete[] matrix[i];

    delete[] matrix;
    delete[] vector;
    delete[] initS;    
    return 0;
}

void SOR(double** matrix, double* vector, double* initS, int n, double omega, double TOL, int maxIter)
{
    int iter;

    double* truesol = new double[n];

    for (int i = 0; i < n; i++)
        truesol[i] = 0;

    for (iter = 0; iter < maxIter; iter++)
    {
        for (int i = 0; i < n; i++)
        {
            double sigma1 = 0.0;
            double sigma2 = 0.0;

            for (int j = 0; j < i; j++)
                sigma1 += matrix[i][j]*truesol[j];

            for (int j = i+1; j < n+1; j++)
                sigma2 += matrix[i][j]*initS[j];

            truesol[i] = (1-omega)*initS[i] + omega*(-sigma1 - sigma2 + vector[i])/matrix[i][i];
        }

        double sumerr = 0.0;

        for (int i = 0; i < n; i++)
        {
            sumerr += pow(truesol[i] - initS[i],2);
        }

        if (sqrt(sumerr) < TOL)
            break;
        
        // Prepare for next iteration
        for (int i = 0; i < n; i++)
            initS[i] = truesol[i];
    }

    cout << "Solution ";
    printVector(truesol, n);
    cout << "Achieved in: " << iter << " iterations." << endl;

    delete[] truesol;
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