#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <SCmathlib.h>

#include <mpi.h>

using namespace std;

int main(int argc, char** argv)

{
    int rank, totalnodes;
  
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int maxT = 10001, nS = 50000;

    double pP, dR, dT, r[3];
    double ProbP, ProbL, ProbE;

    int is, it, ir, n, N;
    int* nt, * Nt;

    pP = 1.0; dR = 1.0; dT = 0.01;

    ProbP = pP*dT;

    nt = new int[maxT];
    Nt = new int[maxT];

    //creating/initialising vector
    for (it = 0;it < maxT;it++) {

        nt[it] = 0; Nt[it] = 0;

    }

    auto start = MPI_Wtime();
    srand(time(NULL) + rank);
  
    int* local_nt = new int[maxT]; int* local_Nt = new int[maxT];
    
    for (int i = 0; i < maxT; i++)
    {
        local_nt[i] = 0;
        local_Nt[i] = 0;
    }

    for (is = 0; is < (nS/totalnodes); is++) 
    {

        N = 0; n = 0;

        for (it = 0;it < (maxT - 1);it++) {

            ProbL = n*dR/dT;

            ProbE = N*(n+1)*dT;

            for (ir = 0; ir < 3; ir++)

                r[ir] = rand() / (double) RAND_MAX;

            if (ProbP >= r[0])

                N = N + 1;

            else if (ProbL >= r[1])

                n = n - 1;

            else if (ProbE >= r[2])

            {
                N = N - 1;
                n = n + 1;
            }

            local_nt[it + 1] = local_nt[it + 1] + n;
            local_Nt[it + 1] = local_Nt[it + 1] + N;
        }
    }
    
    MPI_Allreduce(local_nt, nt, maxT, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
    MPI_Allreduce(local_Nt, Nt, maxT, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    auto end = MPI_Wtime();

    double interval = double(end-start); // time in seconds

    if (rank == 0)
    {
        for (it = 0;it < 20;it++) {

            cout << setiosflags(ios::fixed)

                << setw(8) << setprecision(2)

                << (maxT - (20 - it)) * dT

                << setw(12) << setprecision(5)

                << (double)nt[maxT - (20 - it)] / nS

                << setw(12)

                << (double)Nt[maxT - (20 - it)] / nS

                << endl;
        }

        cout << "Time taken: " << interval << endl;
    }

    MPI_Finalize();

    delete[] nt;
    delete[] Nt;
    delete[] local_Nt;
    delete[] local_nt;

    return 0;

}