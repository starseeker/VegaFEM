ifndef FINITEDIFFERENCETEST
FINITEDIFFERENCETEST=FINITEDIFFERENCETEST

ifndef CLEANFOLDER
CLEANFOLDER=FINITEDIFFERENCETEST
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..

# the object files to be compiled for this utility
FINITEDIFFERENCETEST_OBJECTS=finiteDifferenceTest.o vectorHelper.o

# the libraries this utility depends on
FINITEDIFFERENCETEST_LIBS=elasticForceModel isotropicHyperelasticFEM stvk forceModel sparseMatrix loadList matrix insertRows volumetricMesh getopts graph corotationalLinearFEM clothBW polarDecomposition matrixIO massSpringSystem objMesh configFile minivector imageIO sparseSolver

# the headers in this library
FINITEDIFFERENCETEST_HEADERS=vectorHelper.h

FINITEDIFFERENCETEST_LINK=$(addprefix -l, $(FINITEDIFFERENCETEST_LIBS)) $(SPOOLES_LIB) $(PARDISO_LIB) $(GLUI_LIB) $(IMAGE_LIBS) $(STANDARD_LIBS)

FINITEDIFFERENCETEST_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/finiteDifferenceTest/, $(FINITEDIFFERENCETEST_OBJECTS))
FINITEDIFFERENCETEST_HEADER_FILENAMES=$(addprefix $(R)/utilities/finiteDifferenceTest/, $(FINITEDIFFERENCETEST_HEADERS))
FINITEDIFFERENCETEST_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(FINITEDIFFERENCETEST_LIBS))
FINITEDIFFERENCETEST_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(FINITEDIFFERENCETEST_LIBS))

include $(FINITEDIFFERENCETEST_LIB_MAKEFILES) ../../libraries/glui/Makefile

all: $(R)/utilities/finiteDifferenceTest/finiteDifferenceTest

CURRENT_DIR = $(shell pwd)
$(R)/utilities/finiteDifferenceTest/finiteDifferenceTest: $(FINITEDIFFERENCETEST_OBJECTS_FILENAMES)
	$(CXXLD) $(LDFLAGS) $(FINITEDIFFERENCETEST_OBJECTS) $(FINITEDIFFERENCETEST_LINK) -Wl,-rpath,$(CURRENT_DIR)/$(GLUI_DIR)/lib -o $@; cp $@ $(R)/utilities/bin/

$(FINITEDIFFERENCETEST_OBJECTS_FILENAMES): %.o: %.cpp $(FINITEDIFFERENCETEST_LIB_FILENAMES) $(FINITEDIFFERENCETEST_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(GLUI_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), FINITEDIFFERENCETEST)
clean: cleanfiniteDifferenceTest
endif

deepclean: cleanfiniteDifferenceTest

cleanfiniteDifferenceTest:
	$(RM) $(FINITEDIFFERENCETEST_OBJECTS_FILENAMES) $(R)/utilities/finiteDifferenceTest/finiteDifferenceTest

endif

