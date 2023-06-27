# Gauss_MPI
We parallelize Gauss algorithm using MPI. (The matrix multiplication here is not parallelized and we compute error norm without parallelization. Therefore the entire program is slow but the part corresponding to Gauss algorithm is fast). We obtain speedup by 1.9 for two cores and 3.6 for four cores.
