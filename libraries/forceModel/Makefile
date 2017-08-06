ifndef FORCEMODEL
FORCEMODEL=FORCEMODEL

ifndef CLEANFOLDER
CLEANFOLDER=FORCEMODEL
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this library
FORCEMODELOBJECTS=forceModel.o finiteDifferenceTester.o

# the libraries this library depends on
FORCEMODELLIBS=sparseMatrix

# the headers in this library
FORCEMODELHEADERS=forceModel.h finiteDifferenceTester.h

FORCEMODELOBJECTS_FILENAMES=$(addprefix $(L)/forceModel/, $(FORCEMODELOBJECTS))
FORCEMODELHEADER_FILENAMES=$(addprefix $(L)/forceModel/, $(FORCEMODELHEADERS))
FORCEMODELLIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(FORCEMODELLIBS))
FORCEMODELLIB_FILENAMES=$(call GET_LIB_FILENAMES, $(FORCEMODELLIBS))

include $(FORCEMODELLIB_MAKEFILES)

all: $(L)/forceModel/libforceModel.a

$(L)/forceModel/libforceModel.a: $(FORCEMODELOBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/forceModel/*.h $(L)/include

$(FORCEMODELOBJECTS_FILENAMES): %.o: %.cpp $(FORCEMODELLIB_FILENAMES) $(FORCEMODELHEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), FORCEMODEL)
clean: cleanforceModel
endif

deepclean: cleanforceModel

cleanforceModel:
	$(RM) $(FORCEMODELOBJECTS_FILENAMES) $(L)/forceModel/libforceModel.a

endif

