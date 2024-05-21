#include <iostream> 
#include <iomanip> 
#include <cmath> 
#include <ctime> 
#include <mpi.h>
#include <chrono>

using namespace std; 
using namespace std::chrono;

//const double Pi=4.0*atan(1.0); 

int main(int argc, char* argv[]) 
{
    int points, points_in_circle, points_in_square;
    double approx_pi;

    cout << "Enter number of points: ";
    cin >> points;

    if (points < 100000)
    {
        cout << "Too small! (n >= 100000)";
        return 0;
    }

    srand(time(NULL));

    auto start = high_resolution_clock::now();

    for (int i = 0; i <= points; i++)
    {
        double x = (double) rand()/RAND_MAX;
        double y = (double) rand()/RAND_MAX;

        if ((pow(x, 2) + pow(y, 2)) <= 1) //points INSIDE quarter
        {
            points_in_circle++;
        }
    }

    auto stop = high_resolution_clock::now();

    auto interval = duration<double>(stop-start);

    //cout << points_in_circle << " " << points_in_square << endl;
    approx_pi = (double) 4*(points_in_circle)/points;
    cout << "Approximation of PI value: " << approx_pi << endl;
    cout << "Time taken: " << interval.count() << endl;

} 

 