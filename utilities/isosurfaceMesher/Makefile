ifndef ISOSURFACEMESHER
ISOSURFACEMESHER=ISOSURFACEMESHER

ifndef CLEANFOLDER
CLEANFOLDER=ISOSURFACEMESHER
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for these utilities
ISOSURFACEMESHER_OBJECTS=

# the libraries these utilities depend on
ISOSURFACEMESHER_LIBS=mesher distanceField wildMagic/geometryQuery wildMagic/meshKey wildMagic/rational getopts windingNumber volumetricMesh objMesh graph sparseMatrix imageIO camera matrixIO lighting configFile loadList

# the headers in this library
ISOSURFACEMESHER_HEADERS=

ISOSURFACEMESHER_LINK=$(addprefix -l, $(notdir $(ISOSURFACEMESHER_LIBS))) $(STANDARD_LIBS)

ISOSURFACEMESHER_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/isosurfaceMesher/, $(ISOSURFACEMESHER_OBJECTS))
ISOSURFACEMESHER_HEADER_FILENAMES=$(addprefix $(R)/utilities/isosurfaceMesher/, $(ISOSURFACEMESHER_HEADERS))
ISOSURFACEMESHER_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(ISOSURFACEMESHER_LIBS))
ISOSURFACEMESHER_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(ISOSURFACEMESHER_LIBS))

include $(ISOSURFACEMESHER_LIB_MAKEFILES)

all: $(R)/utilities/isosurfaceMesher/isosurfaceMesher

$(R)/utilities/isosurfaceMesher/isosurfaceMesher: $(R)/utilities/isosurfaceMesher/isosurfaceMesher.cpp $(ISOSURFACEMESHER_LIB_FILENAMES) $(ISOSURFACEMESHER_HEADER_FILENAMES)
	$(CXXLD) $(LDFLAGS) $(INCLUDE) $(GLUI_INCLUDE) $(ISOSURFACEMESHER_OBJECTS) $^ $(ISOSURFACEMESHER_LINK) -Wl,-rpath,$(GLUI_DIR)/lib $(IMAGE_LIBS) -o $@; cp $@ $(R)/utilities/bin/
	
$(ISOSURFACEMESHER_OBJECTS_FILENAMES): %.o: %.cpp $(ISOSURFACEMESHER_LIB_FILENAMES) $(ISOSURFACEMESHER_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(IMAGE_INCLUDE) $(GLUI_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), ISOSURFACEMESHER)
clean: cleanisosurfaceMesher
endif

deepclean: cleanisosurfaceMesher

cleanisosurfaceMesher:
	$(RM) $(R)/utilities/isosurfaceMesher/isosurfaceMesher 

debug:
	make "OPTIMIZATIONFLAG=-g" -j

endif

