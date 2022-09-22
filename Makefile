############################################################
#                         MAKEFILE                         #
############################################################

# GNU Compiler
CMPF  = mpifort -c -cpp -I/home/pietro/petsc/petsc--openmpi/include -pg -g
LNK   = mpifort -cpp -pg -g
#-I/usr/include/suitesparse  -L/usr/lib -lumfpack
#OPTF = -O0 -Wall -Wextra -Warray-temporaries -ggdb3 -pedantic -fimplicit-none -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -mcmodel=medium # Debug options
OPTF = -O3 -Wall -Wextra -fimplicit-none -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -mcmodel=medium # Standard optimization options 
#OPTF = -Ofast -march=native -Wall -Wextra -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -mcmodel=medium # Aggressive optimization options 

# Objects: list of all objects *.o
OBJS = mpi_common.o  global.o  screen.o  tools.o  initialization.o  timecycle.o  grid_and_partition.o  particle.o  collisions.o  postprocess.o  fields.o  umfpack.o  mt19937.o  fully_implicit.o
#OBJDSMC = mpi_common.o tools.o screen.o global.o postprocess.o initialization.o timecycle.o





#  The following variable must either be a path to petsc.pc or just "petsc" if petsc.pc
#  has been installed to a system location or can be found in PKG_CONFIG_PATH.
petsc.pc := $(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig/petsc.pc

# Additional libraries that support pkg-config can be added to the list of PACKAGES below.
PACKAGES := $(petsc.pc)

CC := $(shell pkg-config --variable=ccompiler $(PACKAGES))
CXX := $(shell pkg-config --variable=cxxcompiler $(PACKAGES))
FC := $(shell pkg-config --variable=fcompiler $(PACKAGES))
CFLAGS_OTHER := $(shell pkg-config --cflags-only-other $(PACKAGES))
CFLAGS := $(shell pkg-config --variable=cflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
CXXFLAGS := $(shell pkg-config --variable=cxxflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
FFLAGS := $(shell pkg-config --variable=fflags_extra $(PACKAGES))
CPPFLAGS := $(shell pkg-config --cflags-only-I $(PACKAGES))
LDFLAGS := $(shell pkg-config --libs-only-L --libs-only-other $(PACKAGES))
LDFLAGS += $(patsubst -L%, $(shell pkg-config --variable=ldflag_rpath $(PACKAGES))%, $(shell pkg-config --libs-only-L $(PACKAGES)))
LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm
CUDAC := $(shell pkg-config --variable=cudacompiler $(PACKAGES))
CUDAC_FLAGS := $(shell pkg-config --variable=cudaflags_extra $(PACKAGES))
CUDA_LIB := $(shell pkg-config --variable=cudalib $(PACKAGES))
CUDA_INCLUDE := $(shell pkg-config --variable=cudainclude $(PACKAGES))



# Executable generation by the linker
pantera.exe: pantera.o $(OBJS) 
	$(LNK) $(OPTF) pantera.o $(OBJS) \
	            -o pantera.exe -L/usr/lib -lumfpack $(LDFLAGS) $(LDLIBS) -I/home/pietro/petsc/petsc--openmpi/include

# Objects generation
pantera.o: pantera.f90  $(OBJS) 
	$(CMPF) $(OPTF) pantera.f90

global.o: global.f90  mpi_common.o  particle.o
	$(CMPF) $(OPTF) global.f90

fully_implicit.o: fully_implicit.f90  global.o  screen.o  tools.o
	$(CMPF) $(OPTF) fully_implicit.f90

timecycle.o: timecycle.f90  global.o  particle.o  screen.o  collisions.o  postprocess.o  fields.o  fully_implicit.o
	$(CMPF) $(OPTF) timecycle.f90

initialization.o: initialization.f90  global.o  tools.o  grid_and_partition.o
	$(CMPF) $(OPTF) initialization.f90

tools.o: tools.f90  mpi_common.o  global.o  mt19937.o
	$(CMPF) $(OPTF) tools.f90

grid_and_partition.o: grid_and_partition.f90  mpi_common.o  global.o  tools.o
	$(CMPF) $(OPTF) grid_and_partition.f90

screen.o: screen.f90
	$(CMPF) $(OPTF) screen.f90

mpi_common.o: mpi_common.f90
	$(CMPF) $(OPTF) mpi_common.f90

particle.o: particle.f90
	$(CMPF) $(OPTF) particle.f90

collisions.o: collisions.f90
	$(CMPF) $(OPTF) collisions.f90
	
postprocess.o: postprocess.f90  fields.o
	$(CMPF) $(OPTF) postprocess.f90

fields.o: fields.f90  umfpack.o  fully_implicit.o
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
