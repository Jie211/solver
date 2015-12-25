all:
	nvcc innersolvers.cu main.cu solvers.cu start.cu functions/blas.cu functions/io.cu functions/cudafunc.cu CRS/cg.cu CRS/cr.cu CRS/gcr.cu CRS/kskipcg.cu CRS/kskipcr.cu CRS/vpcg.cu CRS/vpcr.cu CRS/vpgcr.cu -Xcompiler "-fopenmp -O3" -ccbin gcc -use_fast_math -arch=sm_35 -I /usr/local/cuda-7.5/samples/common/inc
