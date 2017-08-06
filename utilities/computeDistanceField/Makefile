ifndef COMPUTEDISTANCEFIELD
COMPUTEDISTANCEFIELD=COMPUTEDISTANCEFIELD

ifndef CLEANFOLDER
CLEANFOLDER=COMPUTEDISTANCEFIELD
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this utility
COMPUTEDISTANCEFIELD_OBJECTS=computeDistanceField.o

# the libraries this utility depends on
COMPUTEDISTANCEFIELD_LIBS=distanceField objMesh graph sparseMatrix matrixIO getopts

# the headers in this library
COMPUTEDISTANCEFIELD_HEADERS=

COMPUTEDISTANCEFIELD_LINK=$(addprefix -l, $(notdir $(COMPUTEDISTANCEFIELD_LIBS))) $(STANDARD_LIBS) #$(SPOOLES_LIB) $(PARDISO_LIB) $(GLUI_LIB) $(IMAGE_LIBS) 

COMPUTEDISTANCEFIELD_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/computeDistanceField/, $(COMPUTEDISTANCEFIELD_OBJECTS))
COMPUTEDISTANCEFIELD_HEADER_FILENAMES=$(addprefix $(R)/utilities/computeDistanceField/, $(COMPUTEDISTANCEFIELD_HEADERS))
COMPUTEDISTANCEFIELD_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(COMPUTEDISTANCEFIELD_LIBS))
COMPUTEDISTANCEFIELD_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(COMPUTEDISTANCEFIELD_LIBS))

include $(COMPUTEDISTANCEFIELD_LIB_MAKEFILES) #../../libraries/glui/Makefile

all: $(R)/utilities/computeDistanceField/computeDistanceField

#CURRENT_DIR = $(shell pwd)
$(R)/utilities/computeDistanceField/computeDistanceField: $(COMPUTEDISTANCEFIELD_OBJECTS_FILENAMES)
	$(CXXLD) $(LDFLAGS) $(COMPUTEDISTANCEFIELD_OBJECTS) $(COMPUTEDISTANCEFIELD_LINK) -o $@; cp $@ $(R)/utilities/bin/

$(COMPUTEDISTANCEFIELD_OBJECTS_FILENAMES): %.o: %.cpp $(COMPUTEDISTANCEFIELD_LIB_FILENAMES) $(COMPUTEDISTANCEFIELD_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(GLUI_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), COMPUTEDISTANCEFIELD)
clean: cleancomputeDistanceField
endif

deepclean: cleancomputeDistanceField

cleancomputeDistanceField:
	$(RM) $(COMPUTEDISTANCEFIELD_OBJECTS_FILENAMES) $(R)/utilities/computeDistanceField/computeDistanceField

endif

