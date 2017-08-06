ifndef TETMESHER
TETMESHER=TETMESHER

ifndef CLEANFOLDER
CLEANFOLDER=TETMESHER
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for these utilities
TETMESHER_OBJECTS=

# the libraries these utilities depend on
TETMESHER_LIBS=mesher wildMagic/geometryQuery wildMagic/meshKey wildMagic/rational getopts windingNumber volumetricMesh objMesh imageIO camera matrixIO lighting configFile loadList

# the headers in this library
TETMESHER_HEADERS=

TETMESHER_LINK=$(addprefix -l, $(notdir $(TETMESHER_LIBS))) $(STANDARD_LIBS)

TETMESHER_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/tetMesher/, $(TETMESHER_OBJECTS))
TETMESHER_HEADER_FILENAMES=$(addprefix $(R)/utilities/tetMesher/, $(TETMESHER_HEADERS))
TETMESHER_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(TETMESHER_LIBS))
TETMESHER_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(TETMESHER_LIBS))

include $(TETMESHER_LIB_MAKEFILES)

all: $(R)/utilities/tetMesher/tetMesher

$(R)/utilities/tetMesher/tetMesher: $(R)/utilities/tetMesher/tetMesher.cpp $(TETMESHER_LIB_FILENAMES) $(TETMESHER_HEADER_FILENAMES)
	$(CXXLD) $(LDFLAGS) $(INCLUDE) $(GLUI_INCLUDE) $(TETMESHER_OBJECTS) $^ $(TETMESHER_LINK) -Wl,-rpath,$(GLUI_DIR)/lib $(IMAGE_LIBS) -o $@; cp $@ $(R)/utilities/bin/
	
$(TETMESHER_OBJECTS_FILENAMES): %.o: %.cpp $(TETMESHER_LIB_FILENAMES) $(TETMESHER_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(IMAGE_INCLUDE) $(GLUI_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), TETMESHER)
clean: cleantetMesher
endif

deepclean: cleantetMesher

cleantetMesher:
	$(RM) $(R)/utilities/tetMesher/tetMesher 

debug:
	make "OPTIMIZATIONFLAG=-g" -j

endif

