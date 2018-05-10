MANSEC           = KSP
CLEANFILES       = rhs.vtk solution.vtk
NP               = 1
CSOURCES 				 = $(wildcard *.cpp)
FCOMPILER        = gfortran -O2

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

f: modules.o porticomb.o $(CSOURCES:.cpp=.o) 
	@-${CLINKER} -o $@ $^ ${PETSC_KSP_LIB} -Wtabs -Wcomment -lboost_system -lgfortran

debug: $(CSOURCES:.cpp=.o)
	@-${CLINKER} -o $@ $^ ${PETSC_KSP_LIB} -g
	@gdb debug

modules.o: modules.for
	@ ${FCOMPILER} -c modules.for

porticomb.o: porticomb.for
	@ ${FCOMPILER} -c porticomb.for

clear:
	@$ rm *.o *~ f *.vtu mirror* domain* *.mod *.dat ma26* tensao* esforc* saida omega.txt

run1:
	@$ mpirun -np 1 ./f

run2:
	@$ mpirun -np 2 ./f

run3:
	@$ mpirun -np 3 ./f

run4:
	@$ mpirun -np 4 ./f

run5:
	@$ mpirun -np 5 ./f

run6:
	@$ mpirun -np 6 ./f

run7:
	@$ mpirun -np 7 ./f

run8:
	@$ mpirun -np 8 ./f

run16:
	@$ mpirun -np 16 ./f -pc_type jacobi -ksp_type gmres -ksp_monitor_singular_value -ksp_gmres_restart 1000
