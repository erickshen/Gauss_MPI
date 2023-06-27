double ErrorNorm(double* x, double* y, int n);
double MatrixNorm(double* x, int n);
double FunctionInput(int k, int n, int i, int j);
int MatrixInputFile(int m, double* x, char* t);
void MatrixOutput(int l, int m, int n, double* x);
void MatrixOutputMPI(double *my_massive,int l,int n,int m, double* secretRow);
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
        double* inverse);
