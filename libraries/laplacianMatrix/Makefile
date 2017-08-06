ifndef LAPLACIANMATRIX
LAPLACIANMATRIX=LAPLACIANMATRIX

ifndef CLEANFOLDER
CLEANFOLDER=LAPLACIANMATRIX
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this library
LAPLACIANMATRIX_OBJECTS=laplacianMatrix.o

# the libraries this library depends on
LAPLACIANMATRIX_LIBS=sparseMatrix volumetricMesh minivector sparseSolver

# the headers in this library
LAPLACIANMATRIX_HEADERS=laplacianMatrix.h

LAPLACIANMATRIX_OBJECTS_FILENAMES=$(addprefix $(L)/laplacianMatrix/, $(LAPLACIANMATRIX_OBJECTS))
LAPLACIANMATRIX_HEADER_FILENAMES=$(addprefix $(L)/laplacianMatrix/, $(LAPLACIANMATRIX_HEADERS))
LAPLACIANMATRIX_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(LAPLACIANMATRIX_LIBS))
LAPLACIANMATRIX_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(LAPLACIANMATRIX_LIBS))

include $(LAPLACIANMATRIX_LIB_MAKEFILES)

all: $(L)/laplacianMatrix/liblaplacianMatrix.a

$(L)/laplacianMatrix/liblaplacianMatrix.a: $(LAPLACIANMATRIX_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/laplacianMatrix/*.h $(L)/include

$(LAPLACIANMATRIX_OBJECTS_FILENAMES): %.o: %.cpp $(LAPLACIANMATRIX_LIB_FILENAMES) $(LAPLACIANMATRIX_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), LAPLACIANMATRIX)
clean: cleanlaplacianMatrix
endif

deepclean: cleanlaplacianMatrix

cleanlaplacianMatrix:
	$(RM) $(LAPLACIANMATRIX_OBJECTS_FILENAMES) $(L)/laplacianMatrix/liblaplacianMatrix.a

endif
