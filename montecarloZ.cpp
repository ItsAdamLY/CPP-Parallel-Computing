#include <iostream> 
#include <iomanip> 
#include <cmath> 
#include <cstdlib> 
#include <ctime> 
#include <random>

using namespace std; 

double fI(double x)
{ 
    return (2/sqrt(M_PI))*pow(M_E, -pow(x, 2)); 
    //return pow(x, 2)*cos(x);
} 


int main() 

{ 

  int i,n,ic=0; 

  double xa,xb,yc,yd,R; 

  double xi,yi,A; 

  srand(time(NULL)); 

  cout << "Enter the number of points, n>=100000: "; 

  cin >> n; 

  xa=0.0; xb=M_PI; yc=0.0; yd=fI(M_PI);//2.0/sqrt(M_PI); 

  R=(xb-xa)*(yd-yc); 

  if (n>=100000) { 

    for(i=1;i<=n;i++) { 

      xi=(double)rand()/(RAND_MAX+1)*(xb-xa)+xa;

      yi=(double)rand()/(RAND_MAX+1)*(yd-yc)+yc;

      if (yi <= fI(xi)) ic++; 

    } 

    A=double(ic)*R/n; 

    cout << setprecision(10); 

    cout << "Integration estimated by " << n << " sample is " << A; 

    cout << endl; 

  } 

  else cout << "Sample size n= " << n << " is so small\n"; 

  return 0; 

} 
