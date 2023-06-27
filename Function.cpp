#include <iostream>
#include <cmath>
double FunctionInput(int k, int n, int i, int j)
{
    if(k == 1)
         return (double) (n - std::max(i+1, j+1) + 1);
    else if(k == 2)
	if(j==1)
	return 0;
	else
        return (double) std::max(i+1, j+1);
    else if(k == 3)
        return (double) sqrt((i-j)*(i-j));
    else if(k == 4)
        return (double) 1/(i+j+1);

    else
        return -1;

}

