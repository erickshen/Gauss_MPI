#include <iostream>
#include <cmath>
double MatrixNorm(double* x, int m)
{
    double Norm = 0;
    double rows = 0;
for(int j = 0; j <= m-1; j++)
{
	Norm = 0;
    for(int i  = 0; i <= m-1; i++)
        Norm += x[j*m+i]*x[j*m+i];
	if(Norm > rows)
		rows = Norm;
}
    return sqrt(rows);
}
