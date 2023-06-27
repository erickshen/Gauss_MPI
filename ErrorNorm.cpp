#include <iostream>
#include <cmath>
double ErrorNorm(double* x, double* y, int n)
{
    std::cout.setf(std::ios::scientific);
    double norm = 0;
    double sum = 0;
        for(int i = 0; i <= n-1; i++)
                for(int k = 0; k <= n-1; k++)
                    {
                        for(int j = 0; j <= n-1; j++)
                            sum+=x[i*n + j%n]*y[j*n + k%n];
                            if(i == k)
                                norm = norm + (sum-1)*(sum-1);
                            else
                                norm = norm + sum*sum;
                            sum = 0;
                    }
        return sqrt(norm);
}
