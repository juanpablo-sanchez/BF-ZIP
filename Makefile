#Makefile to make fast Gibbs sampling code
include ./Makeinit

prog=ZIP_FB_h2_ratios_INV

dir= ./libs

a.out:	$(prog).o  $(dir)/denseop.a $(dir)/stat.a $(dir)/matrix_oper.a
	$(f90)  $(optf90) $(prog).o $(dir)/stat.a $(dir)/matrix_oper.a $(dir)/denseop.a  $(mkl)

	mv a.out ./bin/$(prog)

	       
$(prog).o:	$(prog).f90  
	$(f90)  $(optf90) $(fpp) $(parallel) $(openmp)  $(optdir)$(dir)  -c $(optf90) $(prog).f90

