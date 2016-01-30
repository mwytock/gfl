CFLAGS = -O -largeArrayDims CFLAGS='$$CFLAGS -std=c99'

default:
	mex $(CFLAGS) dtd_factor.c
	mex $(CFLAGS) dtd_solve.c
	mex $(CFLAGS) dtd_inverse.c
	mex $(CFLAGS) tri_factor2.c
	mex $(CFLAGS) tri_solve2.c
	mex $(CFLAGS) tri_inverse2.c
