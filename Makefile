############################################################
#                         MAKEFILE                         #
############################################################

# GNU Compiler
CMP  = mpifort -c -cpp -I$(PETSC_DIR)/include
LNK   = mpifort -cpp
OPTF = -O3 -fimplicit-none

# A few more compiler flags for GNU fortran that could be useful:
#-O0 -Ofast -Wall -Wextra -Warray-temporaries -ggdb3 -pedantic -fimplicit-none -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -mcmodel=medium # Debug options
#-march=native -Wall -Wextra -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -mcmodel=medium # Aggressive optimization options 

# Objects: list of all objects *.o
OBJS = mpi_common.o  global.o  screen.o  tools.o  initialization.o  timecycle.o  grid_and_partition.o  particle.o  collisions.o  postprocess.o  fields.o  mt19937.o  fully_implicit.o


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


all: pantera.exe

debug: CMP += -g
debug: LNK += -g
debug: OPTF += -Wall -Wextra -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow
debug: pantera.exe

# Executable generation by the linker
pantera.exe: pantera.o $(OBJS) 
	$(LNK) $(OPTF) pantera.o $(OBJS) \
	            -o pantera.exe -L/usr/lib $(LDFLAGS) $(LDLIBS) -I$(PETSC_DIR)/include

# Objects generation
pantera.o: pantera.f90  $(OBJS) 
	$(CMP) $(OPTF) pantera.f90

global.o: global.f90  mpi_common.o  particle.o
	$(CMP) $(OPTF) global.f90

fully_implicit.o: fully_implicit.f90  global.o  screen.o  tools.o
	$(CMP) $(OPTF) fully_implicit.f90

timecycle.o: timecycle.f90  global.o  particle.o  screen.o  collisions.o  postprocess.o  fields.o  fully_implicit.o
	$(CMP) $(OPTF) timecycle.f90

initialization.o: initialization.f90  global.o  tools.o  grid_and_partition.o
	$(CMP) $(OPTF) initialization.f90

tools.o: tools.f90  mpi_common.o  global.o  mt19937.o
	$(CMP) $(OPTF) tools.f90

grid_and_partition.o: grid_and_partition.f90  mpi_common.o  global.o  tools.o
	$(CMP) $(OPTF) grid_and_partition.f90

screen.o: screen.f90
	$(CMP) $(OPTF) screen.f90

mpi_common.o: mpi_common.f90
	$(CMP) $(OPTF) mpi_common.f90

particle.o: particle.f90
	$(CMP) $(OPTF) particle.f90

collisions.o: collisions.f90
	$(CMP) $(OPTF) collisions.f90
	
postprocess.o: postprocess.f90  fields.o
	$(CMP) $(OPTF) postprocess.f90

fields.o: fields.f90  fully_implicit.o
	$(CMP) $(OPTF) fields.f90
	
mt19937.o: mt19937.f90
	$(CMP) $(OPTF) mt19937.f90

	
# Cleaning command
clean: 
	@echo cleaning objects, modules and executables 
	rm  -f  *.o  *.mod  *.exe  *~

cleanoutput:
	@echo cleaning output and dump files
	rm  -f  dumps/*

