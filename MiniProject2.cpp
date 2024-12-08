#include <iostream>
#include <mpi.h>
#include <cmath>
#include "JacobiP.h"
#include "Jacobi.h"

using namespace std;

void SOR(double** matrix, double* vector, double* initS, int n, double omega, double TOL, int maxIter);

void Poisson(double AA, int M, int N, double* domain, double* initS, int rank, int nodesize, double TOL, int maxiter);

double f(double A, double x, double y) {return A*sin(2*M_PI*x)*sin(2*M_PI*y);}
double g1(double x) {return 0;}
double g2(double x) {return 0;}
double g3(double y) {return 0;}
double g4(double y) {return 0;}

int main(int argc, char* argv[])
{
    int tprocess, rank;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &tprocess);
    MPI_Comm_rank(comm, &rank);

    double* domain = new double[4];
    double abcd[4] = {0,1,0,1};
    /*     double h = 0.1/4;
    double k = 0.1/4;
    int m = (abcd[1]-abcd[0])/h, n = (abcd[3]-abcd[2])/k; */

    int m = 20, n = 20;

    double* initS = new double[(m+1)*(n+1)];

    for (int i = 0; i < 4; i++)
        domain[i] = abcd[i];

    for (int j = 0; j < (m+1)*(n+1); j++)
        initS[j] = 2;

    auto start = MPI_Wtime();
    Poisson(10, m, n, domain, initS, rank, tprocess, 1e-6, 10000);
    auto end = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);

    double duration = end - start;
    if (rank == 0)
        cout << "Elapsed Time: " << duration << " seconds." << endl;

    delete[] domain;
    delete[] initS;

    MPI_Finalize();
}

void Poisson(double AA, int M, int N, double* domain, double* initS, int rank, int nodesize, double TOL, int maxiter)
{
    double xl = domain[0];
    double xr = domain[1];
    double yb = domain[2];
    double yt = domain[3];

    int m = M + 1, n = N + 1;

    double h = (xr-xl)/M, k = (yt-yb)/N;

    double** A = new double*[m*n];
    double* b = new double[m*n];

    for (int i = 0; i < m*n; i++)
    {
        A[i] = new double[m*n];
        for (int j = 0; j < m*n; j++)
        {
            A[i][j] = 0;
        }
        b[i] = 0;
    }

    double h2 = h*h;
    double k2 = k*k;

    double* x = new double[m];
    double* y = new double[n];

    for (int i = 0; i < m; i++)
        x[i] = xl + i*h;
    
    for (int j = 0; j < n; j++)
        y[j] = yb + j*k;

    // Setting interior points

    for (int i = 1; i < M; i++)
    {
        for (int j = 1; j < N; j++)
        {
            A[i+j*m][i-1+j*m]=1/h2;
            A[i+j*m][i+1+j*m]=1/h2;
            
            A[i+j*m][i+j*m]=-2/h2-2/k2;
            
            A[i+j*m][i+(j-1)*m]=1/k2;
            A[i+j*m][i+(j+1)*m]=1/k2;
            
            b[i+j*m]=f(AA, x[i],y[j]);
        }
    }

    for (int i = 0; i < m; i++)
    {
        A[i][i] = 1;
        b[i]    = g1(x[i]);
        
        A[i+(n-1)*m][i+(n-1)*m] = 1;
        b[i+(n-1)*m] = g2(x[i]);
    }
        
    for (int j = 1; j < N; j++)
    {
        A[j*m][j*m] = 1;
        b[j*m]    = g3(y[j]);
       
        A[m-1+j*m][m-1+j*m] = 1;
        b[m-1+j*m] = g4(y[j]);
    }

    // Solving linear system
    SOR(A, b, initS, m*n, 1.5, 1e-6, 1000000);
    
    for (int i = 0; i < m*n; i++)
        delete[] A[i];
    delete[] A;
    delete[] b;
    delete[] x;
    delete[] y;
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

// To be finished
void SORP(double** matrix, double* vector, double* initS, int n, int rank, int nodesize, double omega, double TOL, int maxIter)
{
    int iter;

    double localN = n/nodesize;

    double* flatMatrix = new double[n*n];
    double* truesol = new double[n];

    // reshape the matrix to be in one vector
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            flatMatrix[i*n + j] = matrix[i][j];
        }
    }

    int remainder = n%nodesize;
    int localNSize;

    if (rank == nodesize - 1 && remainder != 0)
        localNSize = localN + remainder;
    else
        localNSize = localN;

    double* localMat = new double[localNSize*n];
    double* localTrueSol = new double[n/nodesize];

    if (remainder == 0)
    {
        MPI_Scatter(flatMatrix, localN*n, MPI_DOUBLE, localMat, localN*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    else
    {
        int* sendcount = new int[nodesize];
        int* shiftcount = new int[nodesize];

        for (int i = 0; i < nodesize; i++)
        {
            if (i == nodesize - 1)
                sendcount[i] = (localN + remainder)*n;
            
            else
                sendcount[i] = (localN)*n;

            shiftcount[i] = localN*n * i; // by processor
        }

        MPI_Scatterv(flatMatrix, sendcount, shiftcount, MPI_DOUBLE, localMat, localNSize*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        delete[] sendcount;
        delete[] shiftcount;
    }

    double* tempTrueSol = new double[n];

    for (int i = 0; i < n; i++)
        truesol[i] = 0;

    MPI_Bcast(tempTrueSol, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(initS, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (iter = 0; iter < maxIter; iter++)
    {
        for (int i = 0; i < n; i++)
            tempTrueSol[i] = truesol[i];

        // Local matrix vector multiplication 
        for (int i = 0; i < n; i++)
        {
            double* sigma1 = 0.0;
            double sigma2 = 0.0;

            for (int j = 0; j < i; j++)
                temptrue += matrix[i][j]*truesol[j];

            for (int j = i+1; j < n+1; j++)
                sigma2 += matrix[i][j]*initS[j];

            truesol[i] = (1-omega)*initS[i] + omega*(-sigma1 - sigma2 + vector[i])/matrix[i][i];

            delete[] 
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

    delete[] tempTrueSol;
    delete[] localMat;
    delete[] flatMatrix;
    delete[] localTrueSol;
    delete[] truesol;
}