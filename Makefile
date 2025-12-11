#===============================================================================
# 1. Compiler settings
#===============================================================================

#--------------------------------------------------------------
# Fortran compiler
#--------------------------------------------------------------

# Use the MPI parallel HDF5 wrapper for Fortran compiler
FC = h5pfc 

#--------------------------------------------------------------
# Fortran compiler flags
#--------------------------------------------------------------

FFLAGS  = -cpp
FFLAGS += -ffree-line-length-none
FFLAGS += -fdefault-real-8
FFLAGS += -fdefault-double-8
FFLAGS += -fallow-argument-mismatch
FFLAGS += -O3

# Useful for debugging
# FFLAGS += -Og 
# FFLAGS += -g
# FFLAGS += -fbacktrace
# FFLAGS += -Warray-bounds
# FFLAGS += -Wall

#--------------------------------------------------------------
# Advanced Fortran compiler options
#--------------------------------------------------------------

#FFLAGS += -fma -finline-functions
#FFLAGS += -align array64byte
#FFLAGS += -heap-arrays -ip -fno-alias -xHost
#FFLAGS += -xCORE-AVX512 -axCORE-AVX512 -mtune=skylake
#FFLAGS += -axCORE-AVX512
#FFLAGS += -xAVX -axCORE-AVX2

#--------------------------------------------------------------
# System specific libaries
#--------------------------------------------------------------

# Local
# FFTW_LIB = -L/usr/local/lib/

# Common build flags
FFTW_FLAGS = -lfftw3
BLAS_FLAGS = -llapack
HDF5_FLAGS = -lhdf5_fortran -lhdf5 -lz -ldl -lm

LDFLAGS = $(FFTW_FLAGS) $(BLAS_FLAGS) # $(HDF5_FLAGS)

#===============================================================================
# 2. Specify solvers, extensions, modules, boundary and initial conditions
#===============================================================================

#-------------------------------------------------------------------------------
# Object and module directory:
#-------------------------------------------------------------------------------

SRCDIR = src
OBJDIR = obj

#-------------------------------------------------------------------------------
# Target (For production)
#-------------------------------------------------------------------------------

TARGET = afid

#-------------------------------------------------------------------------------
# Boundary conditions, initial conditions and source terms
#-------------------------------------------------------------------------------

# BOUDIR = BoundaryLayer
# BOUDIR = Couette
# BOUDIR = Inflow
# BOUDIR = LidDrivenCavity
BOUDIR = Quiescent

# INIDIR = Couette
# INIDIR = Perturbed
INIDIR = Quiescent
# INIDIR = Uniform

# ADDDIR = ConstFlux
# ADDDIR = ConstShear
ADDDIR = Zero

#-------------------------------------------------------------------------------
# Files that contain modules:
#-------------------------------------------------------------------------------

MFILES = decomp_2d.F90 \
		 Param.F90 \
		 AuxiliaryRoutines.F90 \
		 LagrangianPointParticle.F90

#===============================================================================
# 4. Automatic file discovery and duplicate source files handling
#===============================================================================

SFILES  = $(shell for f in $(SRCDIR)/FlowSolver/*.F90; do [ -f "$$f" ] && basename "$$f"; done)
SFILES += $(shell for f in $(SRCDIR)/InitialConditions/$(INIDIR)/*.F90; do [ -f "$$f" ] && basename "$$f"; done)
SFILES += $(shell for f in $(SRCDIR)/BoundaryConditions/$(BOUDIR)/*.F90; do [ -f "$$f" ] && basename "$$f"; done)
SFILES += $(shell for f in $(SRCDIR)/SourceTerms/$(ADDDIR)/*.F90; do [ -f "$$f" ] && basename "$$f"; done)
SFILES += $(shell for f in $(SRCDIR)/ParticleSolver/*.F90; do [ -f "$$f" ] && basename "$$f"; done)

FFILES = $(shell echo $(SFILES) | tr ' ' '\n' | sort -u)

#-------------------------------------------------------------------------------
# File name substitution
#-------------------------------------------------------------------------------

MOBJS := $(MFILES:%.F90=$(OBJDIR)/%.o)
FOBJS := $(FFILES:%.F90=$(OBJDIR)/%.o)

#===============================================================================
# 5. Make rules
#===============================================================================

all: directories $(TARGET)

$(TARGET): $(MOBJS) $(FOBJS) 
	$(FC) $(FFLAGS) -J $(OBJDIR) -o $@ $^ $(LDFLAGS)

#-------------------------------------------------------------------------------
# Dependencies (Ordering is very important here!)
#-------------------------------------------------------------------------------

# Modules get compiled first

$(OBJDIR)/decomp_2d.o: $(SRCDIR)/2DECOMP/decomp_2d.F90
	$(FC) $(FFLAGS) -J $(OBJDIR) -c $< -o $@ $(LDFLAGS)
$(OBJDIR)/Param.o: $(SRCDIR)/FlowSolver/Param.F90
	$(FC) $(FFLAGS) -J $(OBJDIR) -c $< -o $@ $(LDFLAGS)
$(OBJDIR)/AuxiliaryRoutines.o: $(SRCDIR)/FlowSolver/AuxiliaryRoutines.F90
	$(FC) $(FFLAGS) -J $(OBJDIR) -c $< -o $@ $(LDFLAGS)
$(OBJDIR)/LagrangianPointParticle.o: $(SRCDIR)/ParticleSolver/LagrangianPointParticle.F90
	$(FC) $(FFLAGS) -J $(OBJDIR) -c $< -o $@ $(LDFLAGS)

# Other non-module files get compiled with dependency on module objects

$(OBJDIR)/%.o: $(SRCDIR)/ParticleSolver/%.F90 $(MOBJS)
	$(FC) $(FFLAGS) -J $(OBJDIR) -c $< -o $@ $(LDFLAGS)
$(OBJDIR)/%.o: $(SRCDIR)/InitialConditions/$(INIDIR)/%.F90 $(MOBJS)
	$(FC) $(FFLAGS) -J $(OBJDIR) -c $< -o $@ $(LDFLAGS)
$(OBJDIR)/%.o: $(SRCDIR)/BoundaryConditions/$(BOUDIR)/%.F90 $(MOBJS)
	$(FC) $(FFLAGS) -J $(OBJDIR) -c $< -o $@ $(LDFLAGS)
$(OBJDIR)/%.o: $(SRCDIR)/SourceTerms/$(ADDDIR)/%.F90 $(MOBJS)
	$(FC) $(FFLAGS) -J $(OBJDIR) -c $< -o $@ $(LDFLAGS)
$(OBJDIR)/%.o: $(SRCDIR)/FlowSolver/%.F90 $(MOBJS)
	$(FC) $(FFLAGS) -J $(OBJDIR) -c $< -o $@ $(LDFLAGS)

#-------------------------------------------------------------------------------
# Clean up
#-------------------------------------------------------------------------------
clean:
	rm -rf $(TARGET) $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*genmod* $(OBJDIR)/*.o $(OBJDIR)

.PHONY: directories
directories: $(OBJDIR)
$(OBJDIR):
	mkdir -p ${OBJDIR}