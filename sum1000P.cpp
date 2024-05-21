#include<iostream>
#include<mpi.h>
#include<iomanip>

using namespace std;

int main(int argc, char ** argv)
{
    int mynode, totalnodes, master;
    int startval,endval;
    double sum,accum;

    double precision, start, end;
    MPI_Status status;

    master = 3;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes); // get totalnodes
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode); // get mynode

    precision = MPI_Wtick();
    start = MPI_Wtime();

    sum = 1.0; // zero sum for accumulation

    // Getting start and end values for distribution of tasks wrt processors

    for (int i=mynode;i<=1000;i+=totalnodes)
    {
        sum = sum + i;
    }

    if (mynode!=master)
    {
        MPI_Send(&sum,1,MPI_DOUBLE,master,1,MPI_COMM_WORLD);
    }
    
    else
    {
        for(int j=1;j<totalnodes;j=j+1)
        {
            if (j != master)
            {
                MPI_Recv(&accum,1,MPI_DOUBLE,j,1,MPI_COMM_WORLD, &status);
                sum = sum + accum;
            }
        }
    }
    
    end = MPI_Wtime();

    if(mynode == master)
    {
        cout << setprecision(20);
        cout << "The sum from 1 to 1000 is: " << sum << endl;
    }

    MPI_Finalize();

    return 0;
}