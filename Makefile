############################################################
#                         MAKEFILE                         #
############################################################

# GNU Compiler
CMP  = mpifort -c -cpp -I$(PETSC_DIR)/include
LNK  = mpifort -cpp
OPTF = -O3 -fimplicit-none

# A few more compiler flags for GNU fortran that could be useful:
#-O0 -Ofast -Wall -Wextra -Warray-temporaries -ggdb3 -pedantic -fimplicit-none -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -mcmodel=medium # Debug options
#-march=native -Wall -Wextra -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -mcmodel=medium # Aggressive optimization options 

# Objects: list of all objects *.o
OBJS = $(BUILDDIR)mpi_common.o  $(BUILDDIR)global.o  $(BUILDDIR)screen.o  $(BUILDDIR)tools.o  $(BUILDDIR)initialization.o  $(BUILDDIR)timecycle.o  $(BUILDDIR)grid_and_partition.o  $(BUILDDIR)particle.o  $(BUILDDIR)collisions.o  $(BUILDDIR)postprocess.o  $(BUILDDIR)fields.o  $(BUILDDIR)mt19937.o  $(BUILDDIR)fully_implicit.o  $(BUILDDIR)washboard.o

SRCDIR = src/
BUILDDIR = build/

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
pantera.exe: createbuilddir
pantera.exe: $(BUILDDIR)pantera.o $(OBJS) 
	$(LNK) $(OPTF) $(BUILDDIR)pantera.o $(OBJS) \
	            -o pantera.exe -L/usr/lib $(LDFLAGS) $(LDLIBS) -I$(PETSC_DIR)/include

# Objects generation
$(BUILDDIR)pantera.o: $(SRCDIR)pantera.f90  $(OBJS) createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)pantera.f90

$(BUILDDIR)global.o: $(SRCDIR)global.f90  $(BUILDDIR)mpi_common.o  $(BUILDDIR)particle.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)global.f90

$(BUILDDIR)fully_implicit.o: $(SRCDIR)fully_implicit.f90  $(BUILDDIR)global.o  $(BUILDDIR)screen.o  $(BUILDDIR)tools.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)fully_implicit.f90

$(BUILDDIR)timecycle.o: $(SRCDIR)timecycle.f90  $(BUILDDIR)global.o  $(BUILDDIR)particle.o  $(BUILDDIR)screen.o  $(BUILDDIR)collisions.o  $(BUILDDIR)postprocess.o  $(BUILDDIR)fields.o  $(BUILDDIR)fully_implicit.o  $(BUILDDIR)washboard.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)timecycle.f90

$(BUILDDIR)initialization.o: $(SRCDIR)initialization.f90  $(BUILDDIR)global.o  $(BUILDDIR)tools.o  $(BUILDDIR)grid_and_partition.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)initialization.f90

$(BUILDDIR)tools.o: $(SRCDIR)tools.f90  $(BUILDDIR)mpi_common.o  $(BUILDDIR)global.o  $(BUILDDIR)mt19937.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)tools.f90

$(BUILDDIR)grid_and_partition.o: $(SRCDIR)grid_and_partition.f90  $(BUILDDIR)mpi_common.o  $(BUILDDIR)global.o  $(BUILDDIR)tools.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)grid_and_partition.f90

$(BUILDDIR)screen.o: $(SRCDIR)screen.f90  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)screen.f90

$(BUILDDIR)mpi_common.o: $(SRCDIR)mpi_common.f90  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)mpi_common.f90

$(BUILDDIR)particle.o: $(SRCDIR)particle.f90  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)particle.f90

$(BUILDDIR)collisions.o: $(SRCDIR)collisions.f90  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)collisions.f90
	
$(BUILDDIR)postprocess.o: $(SRCDIR)postprocess.f90  $(BUILDDIR)fields.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)postprocess.f90

$(BUILDDIR)fields.o: $(SRCDIR)fields.f90  $(BUILDDIR)fully_implicit.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)fields.f90
	
$(BUILDDIR)washboard.o: $(SRCDIR)washboard.f90  $(BUILDDIR)tools.o  $(BUILDDIR)global.o  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)washboard.f90

$(BUILDDIR)mt19937.o: $(SRCDIR)mt19937.f90  createbuilddir
	$(CMP) $(OPTF) -o $@ -J$(BUILDDIR) $(SRCDIR)mt19937.f90

	
# Cleaning command
createbuilddir:  clean
	mkdir  build

clean: 
	@echo cleaning objects, modules and executables 
	rm -rf build *.exe

cleanoutput:
	@echo cleaning output and dump files
	rm  -f  dumps/*

