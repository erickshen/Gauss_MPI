#include <iostream>
#include <fstream>
#include <cmath>
int MatrixInputFile(int m, double* x, char* t)
{
    double a;
    std::cout.setf(std::ios::scientific);
    int i = 0;
    std::ifstream Input(t);
    if(!Input)
    {
        std::cout << "Reading error" << std::endl;
        return -1;
    }
    while(Input)
    {

        if(i > m)
            {   std::cout << "Too many arguments " << std::endl;
            return -1;
            }
 
        Input >> x[i];
        if(Input.fail())
        {
		if(Input.eof())
			break;
            Input.clear();
            Input.close();
            std::cout << "Incorrect format" << std::endl;
            return -1;
        }
        i++;
    }


    if(i <= m-1)
        {
            std::cout << "Too few arguments " << std::endl;
            return -1;
        }
	if(i > m)
        {
            std::cout << "Too many arguments " << std::endl;
            return -1;
        }
	//std::cout << i << std::endl;
    Input.close();
    return 0;
}
