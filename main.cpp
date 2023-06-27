#include "mpi.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#include "Definitions.h"
int main(int argc, char* argv[])
{
if(argc < 4)
	{ std::cout << "Wrong number of arguments"; return -1;
}
    int result = 0;
  char* Filename;
	int n = atoi(argv[1]);
	int m = atoi(argv[2]); 
	int k = atoi(argv[3]); //function parameter
	if(m > n || m < 0 || n <= 0)
		{std::cout << "Incorrect arguments" << std::endl; return -1;}
	if(k == 0)
		Filename = argv[4];
    int status = 0;
    int my_rank;
    int total_size;
    int my_row;
    double* massive;
    double* supp;
    double* suppI;
    int* number_of_rows;
    double* inverse;
    int* Column;
    double* Row;
    double* Max;
    int* Index;
    double* secretRowM;
    double* secretRowI;
    int* my_rows_index;
    double time_taken = 0.0;
    double residual = 0.0;
    double elapsed = 0.0;

	MPI_Init(&argc, &argv); //НА СТАРТ
	MPI_Comm_size(MPI_COMM_WORLD, &total_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Barrier(MPI_COMM_WORLD);

	if(total_size >= n && my_rank == 0)
		{
			std::cout << " Слишком много компов... " << std::endl;
			return -1;
		}
	else if(total_size >= n)
		return -1;
		
		

	Row = new double[n];
	Column = new int[n];
	secretRowM = new double[n];
	secretRowI = new double[n];
	
if(my_rank == 0) {
	std::cout << total_size << std::endl;
	Index = new int[total_size];
	Max = new double[total_size];
	number_of_rows = new int[total_size];
	number_of_rows[0] = n/total_size;
	if(n%total_size)
		number_of_rows[0]+=1;
	supp = new double[n*n];
	suppI = new double[n*n];
	 
for(int i = 0; i <= n*n-1; i++)
            {

                if(i/n == i%n)
                    suppI[i] = 1;
                else suppI[i] = 0;
            }
	
if(k == 0)
	{     
	int a = MatrixInputFile(n*n, supp, Filename);
        if(a == -1)
            {std::cout << "Incorrect arguments" << std::endl;
		status = -1;
		MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
		return -1;}
	}	

	
if(k > 0 && k <5) {
        for(int i = 0; i <= n*n-1; i++)
            {
                supp[i] = FunctionInput(k, n, i/n, i%n);
            } }
if(k < 0 || k > 4)
	{
		std::cout << "Wrong argument" << std::endl; 
		status = -1;
		MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
		return -1;
	}

std::cout << "The matrix is of the form: " << std::endl;
	MatrixOutput(m, m, n, supp);
	std::cout << std::endl;

//std::cout << "The inverse matrix is of the form: " << std::endl;
//	MatrixOutput(m, m, n, suppI);
//	std::cout << std::endl;
}
MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD); //End if(my_rank == 0)
if(status == -1)
{
	return -1;
}
//Сначала понять, сколько комó памяти выделить. 

	my_row = 0;
	for(int i = my_rank; i < n; i += total_size)
		my_row++;		//Запишем количество строк
	my_rows_index = new int[my_row];
	for(int i = 0; i < my_row; i++)
		my_rows_index[i] = my_rank + total_size*i; //Запишем действительные номера строк

MPI_Barrier(MPI_COMM_WORLD); 

//Теперь каждомó выделим память и передадим нóлевомó, скока:

massive = new double[my_row*n];
inverse = new double[my_row*n];
MPI_Barrier(MPI_COMM_WORLD); 
MPI_Gather(&my_row, 1, MPI_INT, number_of_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD); 

MPI_Barrier(MPI_COMM_WORLD); 

//Теперь нóлевой всем всё разошлёt
if(my_rank == 0) //себе	
		for(int i = 0; i < my_row; i++)
			for(int j = 0; j < n; j++) {
				massive[i*n +j] = supp[my_rows_index[i]*n+j];
				inverse[i*n+j] = suppI[my_rows_index[i]*n+j];
			}
MPI_Barrier(MPI_COMM_WORLD); 
if(total_size > 1)
for(int i = 0; i < n; i++)
{
	if(i%total_size != 0){
	if(my_rank == 0)
		{
	MPI_Send(supp + i*n, n, MPI_DOUBLE, i%total_size, 0, MPI_COMM_WORLD);
	 }
	else if(my_rank == i%total_size)
		{
	MPI_Recv(massive + (i/total_size)*n, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);}
}
}
MPI_Barrier(MPI_COMM_WORLD); 
//
 MPI_Barrier(MPI_COMM_WORLD); 
//Надеемся, что код выше работает, и ждём завершения

//std::cout << "miracle!" << std::endl;
if(total_size > 1)
for(int i = 0; i < n; i++)
{
	if(i%total_size != 0){
	if(my_rank == 0)
		{//std::cout << " ive sent " << i << " and my rank is " << my_rank << std::endl;
	MPI_Send(suppI + i*n, n, MPI_DOUBLE, i%total_size, 0, MPI_COMM_WORLD);
	 }
	else if(my_rank == i%total_size)
		{ //std::cout << " ive ecieved " << i << "and my rank is " << my_rank << std::endl;
MPI_Recv(inverse + (i/total_size)*n, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);}
}
}
 MPI_Barrier(MPI_COMM_WORLD); 
//Надеемся, что код выше работает, и ждём завершения
//std::cout << "niralce!" << std::endl;
if(my_rank == 0)
	{
		delete[] supp;
		delete[] suppI;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	time_taken = MPI_Wtime();
	int a = GaussMethod(n, total_size, my_rank, my_row, my_rows_index, Max, Index, Column, Row, secretRowM, secretRowI, massive, inverse);
	MPI_Barrier(MPI_COMM_WORLD);
	time_taken = MPI_Wtime() - time_taken;
if(my_rank == 0)
	std::cout << " time taken: " << time_taken << std::endl;
if(a == -1)
{
std::cout << " Process " << my_rank << " is terminated because matrix is singular " << std::endl; return 0;
}
MatrixOutputMPI(inverse, n,n,m, Row);
if(my_rank == 0)
{
   supp = new double[n*n];
   suppI = new double[n*n];

}
	MPI_Barrier(MPI_COMM_WORLD);
	for(int i = 0; i < n; i++)
		{
		if(i%total_size != 0){
			if(my_rank == i%total_size)
				{ //std::cout <<" i send " << i << std::endl;
MPI_Send(inverse + i/total_size*n, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);}
			else if(my_rank == 0) {
				MPI_Recv(suppI + i*n, n, MPI_DOUBLE, i%total_size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //std::cout <<" i recieved " << i << std::endl; 
}
				     }
			
			else if (my_rank == 0)
				{ //std::cout <<"i permute";
				for(int j = 0; j < n; j++)
					suppI[i*n +j] = inverse[i/total_size*n+j];
				}
		}
	
			
	MPI_Barrier(MPI_COMM_WORLD);	

if(my_rank == 0) {
    if(k != 0)
        for(int i = 0; i <= n*n-1; i++)
                supp[i] = FunctionInput(k, n, i/n, i%n);

    if(k == 0)
        MatrixInputFile(n*n, supp, Filename);
    residual = ErrorNorm(suppI, supp, n);
    elapsed = time_taken;
    printf("%s : residual = %e elapsed = %.2f k = %d n = %d m = %d", argv[0], residual, elapsed, k, n, m);
		}
	delete[] massive;
	delete[] inverse;
	delete[] Column;
        delete[] Row;
	delete[] secretRowI;
	delete[] secretRowM;
	if(my_rank == 0)
		{
			delete[] supp;
			delete[] suppI;
		}
	MPI_Finalize();

    return 0;


}

