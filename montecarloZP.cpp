#include <iostream> 
#include <iomanip> 
#include <cmath> 
#include <cstdlib> 
#include <ctime> 
#include <mpi.h>

using namespace std; 

double fI(double x)
{ 
    return (2/M_PI)*pow(M_E, -pow(x, 2)); 
} 

int main(int argc, char* argv[]) 
{ 

  int size, rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int i,n,ic=0; 

  double xa,xb,yc,yd,R; 

  double xi,yi,A; 

  MPI_Status status;

  srand(time(NULL) + rank); 

  if (rank == 0)
  {
    cout << "Enter number of points: ";
    cin >> n;

    for (int i = 1; i <= (size - 1); i++)
    {
        MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  }

  else
    MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

  xa=0.0; xb=1.25; yc=0.0; yd=2.0/sqrt(M_PI); 

  R=(xb-xa)*(yd-yc); 

  int local_ic;

  if (n>=100000) 
  { 

    for(i=1;i<=n/size;i++) { 

      xi=(double)xa + (double)rand() * (xb - xa) / (double)RAND_MAX;

      yi=(double)yc + (double)rand() * (yd - yc) / (double)RAND_MAX;

      if (yi <= fI(xi)) local_ic++; 

    } 

    MPI_Allreduce(&local_ic, &ic, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    A=double(ic)*R/n; 

    if (rank == 0)
    {
        cout << setprecision(10); 

        cout << "Integration estimated by " << n << " sample is " << A; 

        cout << endl; 
    }
  } 

  else 
  {
    cout << "Sample size n= " << n << " is so small\n"; 
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Finalize();

  return 0; 

} 