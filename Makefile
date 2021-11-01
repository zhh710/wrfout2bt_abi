
# Tool-specific substitution variables
FC      = gfortran
FCFLAGS = -mmacosx-version-min=12.0 -fconvert=big-endian  -O3 -fimplicit-none -ffree-form -fno-second-underscore -frecord-marker=4 -funroll-loops -fopenmp -Wall -Wconversion -mieee-fp -fbounds-check -std=f2008  -c
LDFLAGS = -mmacosx-version-min=12.0 -fconvert=big-endian  -O3 -fimplicit-none -ffree-form -fno-second-underscore -frecord-marker=4 -funroll-loops -fopenmp -Wall -Wconversion -mieee-fp -fbounds-check -std=f2008 
LIBS = ${LIBS_CRTM} ${LIBS_NETCDF}
MOD_LIBS = ${MOD_NETCDF} ${MOD_CRTM}
EXE_FILE = wrfout2tb
all:$(EXE_FILE)
include make.filelist
clean:
	-rm *.o *.mod *.a >/dev/null 2>&1
$(EXE_FILE):$(OBJ_FILES)
	$(FC) $(OBJ_FILES) -o $(EXE_FILE) $(LDFLAGS) $(LIBS)
install:
	mv wrfout2tb ./bin
# ...Universal uninstallation
uninstall:
	-rm -fr ./bin >/dev/null 2>&1

# Specify targets that do not generate filesystem objects
.PHONY: all clean install uninstall

# File dependency and compilation rule include files
include make.dependencies
include make.rules
