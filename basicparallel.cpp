#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char** argv) {
    // Initialize MPI environment
    MPI_Init(&argc, &argv);

    // Get the rank of the current process
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Process 0 sends a message to process 1
    if (rank == 0) 
    {
        int message = 42;
        for (int process = 1; process < size; process++)
        {
            //MPI_Send()
            MPI_Send(&message, 1, MPI_INT, process, 0, MPI_COMM_WORLD);
            cout << "Process 0 sent message: " << message << " to Process " << process <<  endl;
        }
    }

    else 
    {
        int message;
        MPI_Recv(&message, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cout << "Process " << rank << " received message: " << message << endl;
        cout << "Address: " << rank << ": " << &message << endl;
    }

    // Finalize MPI environment
    MPI_Finalize();

    return 0;
}
