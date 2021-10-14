############################################################
#                         MAKEFILE                         #
############################################################

# GNU Compiler
CMPF  = mpifort -c -g
LNK   = mpifort
#-I/usr/include/suitesparse  -L/usr/lib -lumfpack
#OPTF = -O0 -Wall -Wextra -Warray-temporaries -ggdb3 -pedantic -fimplicit-none -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -mcmodel=medium # Debug options
#OPTF = -O3 -Wall -Wextra -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow -mcmodel=medium #-fcheck=bounds # Standard optimization options 
OPTF = -Og -Wall -Wextra -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow -mcmodel=medium #-fcheck=bounds # Pietro debug options
#OPTF = -Ofast -march=native -Wall -Wextra -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -mcmodel=medium # Aggressive optimization options 

# Objects: list of all objects *.o
OBJS = mpi_common.o  global.o  screen.o  tools.o  initialization.o  timecycle.o  grid_and_partition.o  particle.o  collisions.o  postprocess.o  fields.o  umfpack.o  mt19937.o
#OBJDSMC = mpi_common.o tools.o screen.o global.o postprocess.o initialization.o timecycle.o

# Executable generation by the linker
pantera.exe: pantera.o $(OBJS) 
	$(LNK) $(OPTF) pantera.o $(OBJS) \
	            -o pantera.exe -L/usr/lib -lumfpack

# Objects generation
pantera.o: pantera.f90 $(OBJS) 
	$(CMPF) $(OPTF) pantera.f90

timecycle.o: timecycle.f90  global.o  particle.o  screen.o  collisions.o  postprocess.o  fields.o
	$(CMPF) $(OPTF) timecycle.f90

initialization.o: initialization.f90  global.o  tools.o  grid_and_partition.o
	$(CMPF) $(OPTF) initialization.f90

tools.o: tools.f90  grid_and_partition.o  mpi_common.o  global.o  mt19937.o
	$(CMPF) $(OPTF) tools.f90

grid_and_partition.o: grid_and_partition.f90 mpi_common.o global.o 
	$(CMPF) $(OPTF) grid_and_partition.f90

global.o: global.f90 mpi_common.o particle.o
	$(CMPF) $(OPTF) global.f90

screen.o: screen.f90
	$(CMPF) $(OPTF) screen.f90

mpi_common.o: mpi_common.f90
	$(CMPF) $(OPTF) mpi_common.f90

particle.o: particle.f90
	$(CMPF) $(OPTF) particle.f90

collisions.o: collisions.f90
	$(CMPF) $(OPTF) collisions.f90
	
postprocess.o: postprocess.f90 fields.o
	$(CMPF) $(OPTF) postprocess.f90

fields.o: fields.f90  umfpack.o
	$(CMPF) $(OPTF) fields.f90

umfpack.o: umfpack.f90
	$(CMPF) $(OPTF) umfpack.f90
	
mt19937.o: mt19937.f90
	$(CMPF) $(OPTF) mt19937.f90
	
# Cleaning command
clean: 
	@echo cleaning objects, modules and executables 
	rm  -f  *.o  *.mod  *.exe  *~

cleanoutput:
	@echo cleaning output and dump files
	rm  -f  dumps/*

# tools.o: tools.f90
# 	$(CMPF) $(OPTF) tools.f90
# 
# screen.o: screen.f90
# 	$(CMPF) $(OPTF) screen.f90
# 
# global.o: global.f90
# 	$(CMPF) $(OPTF) global.f90
# 
# initialization.o: initialization.f90
# 	$(CMPF) $(OPTF) initialization.f90
# 
# timecycle.o: timecycle.f90
# 	$(CMPF) $(OPTF) timecycle.f90
# 
# postprocess.o: postprocess.f90
# 	$(CMPF) $(OPTF) postprocess.f90
# 
# mpi_common.o: mpi_common.f90
# 	$(CMPF) $(OPTF) mpi_common.f90
# 
# # Cleaning command
# clean: 
#	@echo cleaning objects, modules and executables 
#	rm  *.o  *.mod  *.exe  *~
