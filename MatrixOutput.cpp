#include <iostream>
#include <cmath>
#include <stdio.h>
#include "mpi.h"
void MatrixOutput(int l, int m, int n, double* x)
{
	std::cout << std::endl;
    std::cout.setf(std::ios::scientific);
    for(int i = 0; i <= l-1; i++)
    {
        for(int j = 0; j <= m-1; j++)
            std::cout << x[n*i+j] << " " ;
        std::cout << std::endl;
    }
}

void MatrixOutputMPI(double *my_massive,int l,int n,int m, double* secretRow)
{
        int c;
	int my_rank;
	int total_size;
        c = l+m-(abs(l-m)+l+m)/2;
        double* Column;

	        MPI_Comm_size(MPI_COMM_WORLD, &total_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

for(int i = 0; i < c; i++)
        {
                Column = my_massive + (i/total_size)*n;
                for(int j=0; j < n; j++)
                        secretRow[j]=Column[j];
       if(i%total_size != 0)
			{
                        if(my_rank == 0)
                          MPI_Recv(secretRow, n, MPI_DOUBLE, i%total_size,  0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        else if(my_rank == i%total_size)
                                MPI_Send(secretRow, n, MPI_DOUBLE, 0,  0, MPI_COMM_WORLD);
                        }
       
                if(my_rank == 0)
		{
                        for(int j = 0;j < m;j++)
                                printf("%e ",secretRow[j]);
                        printf("\n");
                }
                MPI_Barrier(MPI_COMM_WORLD);
	}
}
