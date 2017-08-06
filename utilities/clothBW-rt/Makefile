ifndef CLOTHBWRT
CLOTHBWRT=CLOTHBWRT

ifndef CLEANFOLDER
CLEANFOLDER=CLOTHBWRT
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this utility
CLOTHBWRT_OBJECTS=clothBW-rt.o initGraphics.o

# the libraries this utility depends on
CLOTHBWRT_LIBS=objMesh imageIO sceneObject objMesh imageIO clothBW integratorSparse integrator sparseSolver elasticForceModel forceModel sparseMatrix insertRows lighting configFile volumetricMesh loadList getopts camera minivector openGLHelper

# the headers in this library
CLOTHBWRT_HEADERS=initGraphics.h

CLOTHBWRT_LINK=$(addprefix -l, $(CLOTHBWRT_LIBS)) $(SPOOLES_LIB) $(BLASLAPACK_LIB) $(GLUI_LIB) $(CG_LIBS) $(IMAGE_LIBS) $(STANDARD_LIBS)

CLOTHBWRT_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/clothBW-rt/, $(CLOTHBWRT_OBJECTS))
CLOTHBWRT_HEADER_FILENAMES=$(addprefix $(R)/utilities/clothBW-rt/, $(CLOTHBWRT_HEADERS))
CLOTHBWRT_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(CLOTHBWRT_LIBS))
CLOTHBWRT_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(CLOTHBWRT_LIBS))

include $(CLOTHBWRT_LIB_MAKEFILES) ../../libraries/glui/Makefile

all: $(R)/utilities/clothBW-rt/clothBW-rt

CURRENT_DIR = $(shell pwd)
$(R)/utilities/clothBW-rt/clothBW-rt: $(CLOTHBWRT_OBJECTS_FILENAMES)
	$(CXXLD) $(LDFLAGS) $(CLOTHBWRT_OBJECTS) $(CLOTHBWRT_LINK) -Wl,-rpath,$(CURRENT_DIR)/$(GLUI_DIR)/lib -o $@; cp $@ $(R)/utilities/bin/

$(CLOTHBWRT_OBJECTS_FILENAMES): %.o: %.cpp $(CLOTHBWRT_LIB_FILENAMES) $(CLOTHBWRT_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(CG_INCLUDE) $(GLUI_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), CLOTHBWRT)
clean: cleanclothBW-rt
endif

deepclean: cleanclothBW-rt

cleanclothBW-rt:
	$(RM) $(CLOTHBWRT_OBJECTS_FILENAMES) $(R)/utilities/clothBW-rt/clothBW-rt

endif

