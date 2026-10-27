FC = gfortran

FDM_FLAGS = -O2 -g -fbacktrace -ffree-line-length-none
MMM_FLAGS = -O2 -g -fbacktrace -ffree-line-length-none

# ============================================================
# Default target
# ============================================================

all: FDM_1Dconvect.exe MMM_1Dconvect.exe

# ============================================================
# FDM solver
# ============================================================

FDM_1Dconvect.exe: FDM_module.f90 FDM_1Dconvect.f90
	$(FC) $(FDM_FLAGS) -o $@ FDM_module.f90 FDM_1Dconvect.f90

# ============================================================
# MMM solver
# ============================================================

MMM_1Dconvect.exe: MMM_module.f90 MMM_module_TimeINT.f90 MMM_1Dconvect.f90
	$(FC) $(MMM_FLAGS) -o $@ \
	    MMM_module.f90 \
	    MMM_module_TimeINT.f90 \
	    MMM_1Dconvect.f90

# ============================================================
# Clean
# ============================================================

.PHONY: all clean

clean:
	rm -f *.exe *.mod *.o


