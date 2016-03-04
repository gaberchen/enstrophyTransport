# enstrophyTransport
OpenFOAM utility to compute the terms of the enstrophy tranport equation

To compile using wmake, ensure the following files are in the Make folder:

files

options

The "files" file must contain the following:

enstrophyFoam.C

EXE = $(FOAM_USER_APPBIN)/enstrophyFoam

The "options" file must contain the following:

EXE_INC = \
    -I$(LIB_SRC)/postProcessing/postCalc \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude

EXE_LIBS = \
    $(FOAM_LIBBIN)/postCalc.o \
    -lgenericPatchFields \
    -lfiniteVolume \
    -lmeshTools

