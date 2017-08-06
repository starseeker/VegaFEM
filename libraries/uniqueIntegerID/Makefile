ifndef UNIQUEINTEGERID
UNIQUEINTEGERID=UNIQUEINTEGERID

ifndef CLEANFOLDER
CLEANFOLDER=UNIQUEINTEGERID
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header


# the object files to be compiled for this library
UNIQUEINTEGERID_OBJECTS=uniqueIntegerID.o

# the libraries this library depends on
UNIQUEINTEGERID_LIBS=

# the headers in this library
UNIQUEINTEGERID_HEADERS=uniqueIntegerID.h


UNIQUEINTEGERID_OBJECTS_FILENAMES=$(addprefix $(L)/uniqueIntegerID/, $(UNIQUEINTEGERID_OBJECTS))
UNIQUEINTEGERID_HEADER_FILENAMES=$(addprefix $(L)/uniqueIntegerID/, $(UNIQUEINTEGERID_HEADERS))
UNIQUEINTEGERID_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(UNIQUEINTEGERID_LIBS))
UNIQUEINTEGERID_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(UNIQUEINTEGERID_LIBS))

include $(UNIQUEINTEGERID_LIB_MAKEFILES)

all: $(L)/uniqueIntegerID/libuniqueIntegerID.a

$(L)/uniqueIntegerID/libuniqueIntegerID.a: $(UNIQUEINTEGERID_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/uniqueIntegerID/*.h $(L)/include

$(UNIQUEINTEGERID_OBJECTS_FILENAMES): %.o: %.cpp $(UNIQUEINTEGERID_LIB_FILENAMES) $(UNIQUEINTEGERID_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), UNIQUEINTEGERID)
clean: cleanuniqueIntegerID
endif

deepclean: cleanuniqueIntegerID

cleanuniqueIntegerID:
	$(RM) $(UNIQUEINTEGERID_OBJECTS_FILENAMES) $(L)/uniqueIntegerID/libuniqueIntegerID.a

endif
