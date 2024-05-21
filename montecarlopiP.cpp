#include <iostream> 
#include <iomanip> 
#include <cmath> 
#include <stdlib.h> 
#include <ctime> 
#include <mpi.h>

using namespace std; 

//const double Pi=4.0*atan(1.0); 

int main(int argc, char* argv[]) 
{
    int size, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int points, points_in_circle = 0;
    double approx_pi, x, y;

    MPI_Status status;

    if (rank == 0)
    {
        cout << "Enter number of points: ";
        cin >> points;

        if (points < 100000)
        {
            cout << "Too small! (n >= 100000)" << endl;;
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 0;
        }

        for (int i = 1; i <= (size - 1); i++)
        {
            MPI_Send(&points, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }

    else
    {
        MPI_Recv(&points, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }

    double start = MPI_Wtime();

    srand(time(NULL) + rank); // to ensure numbers generate differently for each processes

    int local_points_in_circle = 0;

    for (int i = 0; i < (points/size); i++) //distribute points into several processors
    {
        double x = (double) rand()/RAND_MAX;
        double y = (double) rand()/RAND_MAX;

        if ((pow(x,2) + pow(y,2)) <= 1)
          local_points_in_circle++;
    }

    /*if (rank != 0)
      MPI_Send(&local_points_in_circle, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    
    else
    {
        int temp_points_in_circle;
        points_in_circle = local_points_in_circle;
        for (int i = 1; i < size; i++)
        {
            MPI_Recv(&temp_points_in_circle, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            points_in_circle += temp_points_in_circle;
        }
    }*/

    MPI_Allreduce(&local_points_in_circle, &points_in_circle, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    double stop = MPI_Wtime();
    
    if (rank == 0)
    {
        approx_pi = (double) 4*(points_in_circle)/points;
        cout << "Approximation of PI value: " << approx_pi << endl;
        setprecision(20);
        cout << "Time taken: " << stop-start << endl;
    }

    MPI_Finalize();

    return 0;
} 