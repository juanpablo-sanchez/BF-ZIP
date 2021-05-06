# Makes libraries

include ../Makeinit

all: denseop.a matrix_oper.a stat.a

denseop.a: lapack90r.o  kind.o denseop.o 
	ar cr  denseop.a kind.o denseop.o lapack90r.o

stat.a: ranlib.o prob.o 
	ar cr stat.a prob.o  ranlib.o

matrix_oper.a: matrix_oper.o matrix_oper_sub.o 
	ar cr matrix_oper.a matrix_oper.o matrix_oper_sub.o 
	
kind.o:	kind.f90
	$(f90) -c  $(optf90) kind.f90

denseop.o:	denseop.f90
	$(f90) -c  $(optf90)  denseop.f90

lapack90r.o:	lapack90r.f90 
	$(f90) -c  $(optf90) lapack90r.f90

prob.o:	prob.f90
	$(f90) -c  $(optf90) prob.f90
	
ranlib.o:	ranlib.f
	$(f90) -c  $(optf77) ranlib.f

matrix_oper.o:	matrix_oper.f90 
	$(f90) -c  $(optf90) matrix_oper.f90 

matrix_oper_sub.o:   matrix_oper_sub.F90
	$(f90) -c  $(optf90) $(fpp) $(parallel) $(openmp) -D _MYMKL=1 matrix_oper_sub.F90 

clean:
	rm *.o *.a *.mod


