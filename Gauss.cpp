#include <iostream>
#include <cmath>
#include <fstream>
#include "Definitions.h"
#include <pthread.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include "sync.h"
#include "mpi.h"

int Rationing(double* x, double* inverse, int n);
void PermuteColumn(double* x, int index, int step, double*Row, int n);
int FindMax(double* x, int step, int n);
void PermuteRow(double* x, int index, int step, double*Row, int n);
void Substitute(double* x, double* inverse, int step, int n);
int Back(double* x, double* inverse, int n);
int Rationing(double* x, double* inverse, int n);
double mach_eps(void);
static double eps;
int GaussMethod(	
	int n,
	int total_size,
	int my_rank,
	double my_rows,
	int* my_rows_index,
	double* maxx,
	int* indexx,
	int* Column,
        double* Row, //Для перестановки столбцов
	double* secretRowM,
	double* secretRowI,
        double* my_massive,
        double* inverse)
{
	double max = 0;
	double true_max = 0;
	
	int index;
	int status = 0;

	 eps = 10*mach_eps();
	MPI_Barrier(MPI_COMM_WORLD);



//Нормировка.
	for(int i = 0; i < my_rows; i++)
		{
			double Norm = 0;
			for(int j = 0; j <n; j++)
				Norm += sqrt(my_massive[i*n+j]*my_massive[i*n+j]);
			if(Norm > max)
				max = Norm;
		}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&max, &true_max, 1, MPI_DOUBLE,  MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&true_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	for(int i = 0; i < my_rows; i++)
		for(int j = 0; j < n; j++)
			{
				my_massive[i*n+j] /= true_max;
				inverse[i*n+j] /=true_max;
			}







	MPI_Barrier(MPI_COMM_WORLD);
	int prestep = 0;
	int invstep = 0;
        for(int step = 0; step <= n-1; step++)
    { //Поиск максимóма
 	MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0)
			for(int i = 0; i<total_size; i++)
				maxx[i] = 0.0;
	max = 0.0;
	MPI_Barrier(MPI_COMM_WORLD);
    	for(int i = 0; i < my_rows; i++)
      		{ 
		if(my_rows_index[i] >= step) 
		  for(int j = step; j < n; j++)
	  		{
			if(sqrt(my_massive[i*n + j]*my_massive[i*n + j]) > sqrt(max*max))
        		{
           			max = my_massive[i*n+j];
            		 	index = my_rows_index[i]*n+j;
        		}	
						
      		} 
	} 
//		//Пересылаем нóлевомó максимóмы и индексы
	MPI_Barrier(MPI_COMM_WORLD); 
		MPI_Gather(&max, 1, MPI_DOUBLE, maxx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD); 
		MPI_Gather(&index, 1, MPI_INT, indexx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD); 
		//Всё переслалося... Выбираем максимóм
	if(my_rank == 0)
       	{
		//for(int p =0; p <total_size;p++)
		//	std::cout << " " << maxx[p] << std::endl;
		int ind_max = n*n;
		double k = 0.0;
		
		for(int i = 0; i<total_size; i++)
			{
			if(sqrt(maxx[i]*maxx[i]) > sqrt(k*k))
				{ 
				   k = maxx[i];
				   ind_max = indexx[i]; 
				}
			else if(fabs(sqrt(maxx[i]*maxx[i]) - sqrt(k*k))<eps/100 && indexx[i] < ind_max)
				ind_max = indexx[i];	
			}
		
		//Выбрали максимóм. Теперь надо как-то переставить
		Column[step] = ind_max;
		if(sqrt(k*k) < eps)
			status = -1;
	//std::cout << "\n maximum at the step " << step << " is " <<k << "\n" << std::endl;
	//std::cout << "\n index is " << ind_max << "\n" << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD); 		
		MPI_Bcast(&Column[step], 1, MPI_INT, 0, MPI_COMM_WORLD); //Теперь все знают имя этого героя!
		MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		MPI_Barrier(MPI_COMM_WORLD);
		if(status == -1)
			{
			//std::cout << " process " << my_rank << " terminating " << std::endl;
			 return -1; 
			}


		if(my_rank == step%total_size)
			{
				for(int j = 0; j < n; j++)
				{	secretRowM[j] = my_massive[step/total_size*n + j];
					secretRowI[j] = inverse[step/total_size*n+j];
				}
			}
		else if (my_rank == (Column[step]/n)%total_size)
			{
				for(int j = 0; j < n; j++)
				{
					secretRowM[j] = my_massive[(Column[step]/n)/total_size*n + j];
					secretRowI[j] = inverse[(Column[step]/n)/total_size*n + j];
				}
			}
		MPI_Barrier(MPI_COMM_WORLD);
if((Column[step]/n)%total_size != step%total_size) // Если они разные, то пересылка
{ 
		if(my_rank == (Column[step]/n)%total_size)
			MPI_Send(secretRowM, n, MPI_DOUBLE, step%total_size, 0, MPI_COMM_WORLD);
		else if(my_rank == step%total_size)
			MPI_Recv(secretRowM, n, MPI_DOUBLE, (Column[step]/n)%total_size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		MPI_Barrier(MPI_COMM_WORLD);

		if(my_rank == step%total_size)
			MPI_Send(my_massive + step/total_size*n, n, MPI_DOUBLE, (Column[step]/n)%total_size, 0, MPI_COMM_WORLD);
		else if(my_rank == (Column[step]/n)%total_size)
			MPI_Recv(secretRowM, n, MPI_DOUBLE, step%total_size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		MPI_Barrier(MPI_COMM_WORLD);
			if(my_rank == step%total_size)
				for(int j = 0; j < n; j++)
					my_massive[step/total_size*n + j] = secretRowM[j];
			else if(my_rank == (Column[step]/n)%total_size)
				for(int j = 0; j < n; j++)
					my_massive[(Column[step]/n)/total_size*n + j] = secretRowM[j];
}

else //иначе просто переставим
{
if(my_rank == step%total_size)
		for(int j = 0; j < n; j++)
			{ 
		secretRowM[j] = my_massive[(Column[step]/n)/total_size*n +j];
		my_massive[(Column[step]/n)/total_size*n +j] = my_massive[step/total_size*n+j];
		my_massive[step/total_size*n+j] = secretRowM[j];
			}
		

}// MPI_Barrier(MPI_COMM_WORLD);  MPI_Barrier(MPI_COMM_WORLD);
//Пересылка массивов кончилась, дальше пересылка обратных массивов:

		MPI_Barrier(MPI_COMM_WORLD);
if((Column[step]/n)%total_size != step%total_size) //Если они разные, то пересылка
{
		if(my_rank == (Column[step]/n)%total_size)
			MPI_Send(secretRowI, n, MPI_DOUBLE, step%total_size, 0, MPI_COMM_WORLD);
		else if(my_rank == step%total_size)
			MPI_Recv(secretRowI, n, MPI_DOUBLE, (Column[step]/n)%total_size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		MPI_Barrier(MPI_COMM_WORLD);

		if(my_rank == step%total_size)
			MPI_Send(inverse + step/total_size*n, n, MPI_DOUBLE, (Column[step]/n)%total_size, 0, MPI_COMM_WORLD);
		else if(my_rank == (Column[step]/n)%total_size)
			MPI_Recv(secretRowI, n, MPI_DOUBLE, step%total_size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


		MPI_Barrier(MPI_COMM_WORLD);
			if(my_rank == step%total_size)
				for(int j = 0; j < n; j++)
					inverse[step/total_size*n + j] = secretRowI[j];
			else if(my_rank == (Column[step]/n)%total_size)
				for(int j = 0; j < n; j++)
					inverse[(Column[step]/n)/total_size*n + j] = secretRowI[j];
}

else //иначе просто переставим
{	if(my_rank == step%total_size)
		for(int j = 0; j < n; j++)
			{
		secretRowI[j] = inverse[(Column[step]/n)/total_size*n +j];
		inverse[(Column[step]/n)/total_size*n +j] = inverse[step/total_size*n+j];
		inverse[step/total_size*n+j] = secretRowI[j];
			}
		
}  
//Переставили строчки			

		MPI_Barrier(MPI_COMM_WORLD);
		for(int i = 0; i < my_rows; i++) //Переставим в КАЖДОМ столбцы
		 	{				
			Row[i] = my_massive[i*n + Column[step]%n];
			my_massive[i*n + Column[step]%n] = my_massive[i*n + step%n];
			my_massive[i*n+step%n] = Row[i];
			}
		MPI_Barrier(MPI_COMM_WORLD);
	 //sdkjfsiufsduhfsdjfsjdfsdhsdfjhsdfkjhsfkjsfkjskjfh
	//
	//Рассылка строчек перед вычитанием
	MPI_Barrier(MPI_COMM_WORLD);
	if(step == my_rows_index[prestep])
	{
		MPI_Bcast(&my_massive[prestep*n], n, MPI_DOUBLE, my_rank, MPI_COMM_WORLD);
		for(int i = 0; i < n; i++)
			secretRowM[i] = my_massive[prestep*n+i];
		prestep++;
	}
	else
	 MPI_Bcast(secretRowM, n, MPI_DOUBLE, step%total_size, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if(step == my_rows_index[prestep])
	{
		MPI_Bcast(&inverse[prestep*n], n, MPI_DOUBLE, my_rank, MPI_COMM_WORLD);
		for(int i = 0; i < n; i++)
			secretRowI[i] = inverse[prestep*n+i];
		prestep++;
	}
	else
	 MPI_Bcast(secretRowI, n, MPI_DOUBLE, step%total_size, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	//Вычитание строчек

	for(int j = prestep; j < my_rows; j++)
	{ double a = my_massive[j*n+step]/secretRowM[step];
		for(int i = step; i < n; i++)
			my_massive[j*n+i] = my_massive[j*n+i] - secretRowM[i]*a;
		for(int i = 0; i < n; i++)
			inverse[j*n+i] = inverse[j*n+i] - secretRowI[i]*a;
	}

	
 MPI_Barrier(MPI_COMM_WORLD); 

	
    }



























	MPI_Barrier(MPI_COMM_WORLD);
	//std::cout << std::endl;
	//sleep(2);
	//
	//sleep(2);

MPI_Barrier(MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
for(int i = n-1; i>= 0; i--)
{
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank == i%total_size)
		{for(int j =0; j<n;j++)
			inverse[(i/total_size)*n + j] /= my_massive[(i/total_size)*n +i];
		for(int j = 0; j <n; j++)
			secretRowI[j] = inverse[(i/total_size)*n + j];
		MPI_Bcast(inverse+(i/total_size)*n, n, MPI_DOUBLE, my_rank, MPI_COMM_WORLD);}
	else
		MPI_Bcast(secretRowI, n, MPI_DOUBLE, i%total_size, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	for(int j = 0; j < my_rows; j++)
		{
			if(my_rank != i%total_size || i/total_size != j)
				for(int k = 0; k <n; k++)
					inverse[j*n +k] -= secretRowI[k]*my_massive[j*n+i];
		}
	MPI_Barrier(MPI_COMM_WORLD);
	
MPI_Barrier(MPI_COMM_WORLD); //std::cout<<std::endl;
}
/*MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank ==  && i == n-1) {sleep(2);
				for(int j = 0; j < n; j++)
					std::cout << " " << Column[j] << " " << std::endl;
		sleep(2);} */
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	//Осталось переставить строчки в обратной матрице.
	//std::cout << "am here" << std::endl;
/*
std::cout << std::endl;
sleep(2);

std::cout << std::endl; */
  for(int i = n-1; i >= 0; i--)
	{//if(my_rank == 0) std::cout << "permuting rows : " << i << " and " << (Column[i]%n) << std::endl;
MPI_Barrier(MPI_COMM_WORLD);
	/*if(my_rank == 1 && i == n-1) {sleep(2);
				for(int j = 0; j < n; j++)
					std::cout << " " << Column[j] << " " << std::endl;
		sleep(2);} */ 

		//Переставить строчки с индексами Column[i]%n и i
if(i%total_size != (Column[i]%n)%total_size)
{	//std::cout << "am in" << std::endl;
		if(my_rank == i%total_size)
			{
			MPI_Send(inverse + (i/total_size)*n, n, MPI_DOUBLE, (Column[i]%n)%total_size, 0, MPI_COMM_WORLD);
			MPI_Recv(secretRowI, n, MPI_DOUBLE, (Column[i]%n)%total_size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for(int j = 0; j < n; j++)
				inverse[(i/total_size)*n + j] = secretRowI[j];
			}
		
		else if(my_rank == (Column[i]%n)%total_size)
			{ 
			MPI_Recv(secretRowI, n, MPI_DOUBLE, i%total_size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(inverse+((Column[i]%n)/total_size)*n, n, MPI_DOUBLE, i%total_size, 0, MPI_COMM_WORLD);
				for(int j = 0; j < n; j++)
					inverse[((Column[i]%n)/total_size)*n + j] = secretRowI[j];
			
			}	
}
else if (my_rank==i%total_size)
    for(int j = 0; j < n; j++)
		{
      			secretRowI[j] = inverse[(i/total_size)*n +j];
			inverse[i/total_size*n +j] = inverse[((Column[i]%n)/total_size)*n+j];
			inverse[((Column[i]%n)/total_size)*n+j] = secretRowI[j];
		}
/*std::cout << std::endl;
sleep(2);

std::cout << std::endl;		
*/		
	}
	 
	return 0;

}

double mach_eps(void)
{ double eps = 1.0;
 while(1.0 + eps/2.0 > 1.0)
{ eps /= 2.0;}
return eps;
}
