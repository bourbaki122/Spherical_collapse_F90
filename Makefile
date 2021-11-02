FF90 = gfortran
CCCC = gcc
USR=/usr
LIBS_fortran = -lgfortran -llapack -lblas 
LIBS_fgsl = -L/usr/local/lib -lfgsl -lgsl -lgslcblas -lm
WARNINGS = -fcheck=all -Wall -Wextra -Wconversion

#PREDEFINED VARIABLES
build:
	$(FF90) -c -o rkf45.o rkf45.f90
	$(FF90) -c -o zoomin.o zoomin.f90
	$(FF90) -c -o nms.o nms.f90
	$(FF90) -c -o spline.o spline.f90


comp:
	$(FF90) -g $(WARNINGS) -o SC_IMPLEMENTATION.out -I/usr/local/include/fgsl C_models.f90 SC_IMPLEMENTATION.f90 rkf45.o zoomin.o nms.o spline.o $(LIBS_fgsl)

run:
	./SC_IMPLEMENTATION.out

clean:
	rm -rf  *.mod *.eps *.out
