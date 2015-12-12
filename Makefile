all:
	nvcc innersolvers.c main.c solvers.c start.c functions/blas.c functions/io.c functions/cudafunc.cu CRS/cg.c CRS/cr.c CRS/gcr.c CRS/kskipcg.c CRS/kskipcr.c CRS/vpcg.c CRS/vpcr.c CRS/vpgcr.c -Xcompiler "-fopenmp " -ccbin gcc
