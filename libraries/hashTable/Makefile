ifndef INSERTROWS
INSERTROWS=INSERTROWS

ifndef CLEANFOLDER
CLEANFOLDER=INSERTROWS
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..


# the object files to be compiled for this library
INSERTROWS_OBJECTS=hashTable.o

# the libraries this library depends on
INSERTROWS_LIBS=

# the headers in this library
INSERTROWS_HEADERS=hashTable.h


INSERTROWS_OBJECTS_FILENAMES=$(addprefix $(L)/hashTable/, $(INSERTROWS_OBJECTS))
INSERTROWS_HEADER_FILENAMES=$(addprefix $(L)/hashTable/, $(INSERTROWS_HEADERS))
INSERTROWS_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(INSERTROWS_LIBS))
INSERTROWS_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(INSERTROWS_LIBS))

include $(INSERTROWS_LIB_MAKEFILES)

all: $(L)/hashTable/libhashTable.a

$(L)/hashTable/libhashTable.a: $(INSERTROWS_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/hashTable/*.h $(L)/include

$(INSERTROWS_OBJECTS_FILENAMES): %.o: %.cpp $(INSERTROWS_LIB_FILENAMES) $(INSERTROWS_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), INSERTROWS)
clean: cleanhashTable
endif

deepclean: cleanhashTable

cleanhashTable:
	$(RM) $(INSERTROWS_OBJECTS_FILENAMES) $(L)/hashTable/libhashTable.a

endif
