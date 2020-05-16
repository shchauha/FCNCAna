#
# stuff to make
#
ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else

all: build
help:
	@echo "Available Targets:";\
	cat Makefile | perl -ne 'printf("\t%-15s %s\n",$$1,$$2) if(/^(\S+):[^#]+(#.*)$$/)'

ifndef VERBOSE
  QUIET := @
endif

opt: CFLAGS+=-O2
small: CFLAGS+=-Os
debug: CFLAGS+=-g -O0

#ROOFITINCLUDE = 
#ifdef CMSSW_VERSION
#	ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
#endif


CC = g++
CMSROOT = ./
INCLUDE = $(shell root-config --cflags) -I$(CMSROOT) -I$(CMSROOT)/CORE
# CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(INCLUDE) $(EXTRACFLAGS)
CFLAGS =  -Wall -Wno-unused-function -fPIC -Wno-sign-compare -Wno-unused-variable $(INCLUDE) $(EXTRACFLAGS)
# CFLAGS =  -Wall -Wno-unused-function -g -O1 -fPIC -Wno-sign-compare -Wno-unused-variable $(INCLUDE) $(EXTRACFLAGS)
# CFLAGS =  -Wall -g -O0 -fPIC $(INCLUDE) $(EXTRACFLAGS) -Wall -Wextra
# CFLAGS = -Wall -Wno-unused-function -g -O0 -fPIC  $(INCLUDE) $(EXTRACFLAGS) # default
# CFLAGS = -Wall -Wno-unused-function -g -O1 -fPIC $(INCLUDE) $(EXTRACFLAGS)
ROOTLIBS = $(shell root-config --ldflags --cflags --libs) -lTMVA #-lEG -lGenVector
COREDIR = CORE

# SimpleJetCorrector.icc is horribly coded, so need to force it to be optimized
$(COREDIR)/Tools/JetCorrector.o : CFLAGS+=-O3

DICTINCLUDE = $(ROOTSYS)/include/Math/QuantFuncMathCore.h $(ROOTSYS)/include/TLorentzVector.h $(ROOTSYS)/include/Math/Vector4D.h

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags --libs) -lEG -lGenVector -lTMVA

#DIR = /Users/cerati/SSAnalysis/SSAnalysis/
DIR = ./

CORESOURCES=$(DIR)/$(COREDIR)/CMS3.cc \
 $(DIR)/$(COREDIR)/Base.cc \
 $(DIR)/$(COREDIR)/Config.cc \
 $(DIR)/$(COREDIR)/ElectronSelections.cc \
 $(DIR)/$(COREDIR)/MuonSelections.cc \
 $(DIR)/$(COREDIR)/JetSelections.cc \
 $(DIR)/$(COREDIR)/MetSelections.cc \
 $(DIR)/$(COREDIR)/VertexSelections.cc \
 $(DIR)/$(COREDIR)/TriggerSelections.cc \
 $(DIR)/$(COREDIR)/MCSelections.cc \
 $(DIR)/$(COREDIR)/IsolationTools.cc \
 $(DIR)/$(COREDIR)/Tools/MT2/MT2Utility.cc \
 $(DIR)/$(COREDIR)/Tools/MT2/MT2.cc \
 $(DIR)/$(COREDIR)/Tools/utils.cc \
 $(DIR)/$(COREDIR)/Tools/JetResolution/JetResolution.cxx \
 $(DIR)/$(COREDIR)/Tools/JetCorrector.cc \
 $(DIR)/$(COREDIR)/Tools/goodrun.cc \
 $(DIR)/$(COREDIR)/Tools/datasetinfo/getDatasetInfo.cc \
 $(DIR)/$(COREDIR)/Tools/btagsf/BTagCalibrationStandalone.cc \
 $(DIR)/$(COREDIR)/SSSelections.cc

COREOBJECTSTMP=$(CORESOURCES:.cc=.o) # replace .cc with .o
COREOBJECTS=$(COREOBJECTSTMP:.cxx=.o) # replace .cxx with .o
# $(info COREOBJECTS="$(COREOBJECTS)") # print out
CORELIB=libCMS3CORE.so

SOURCES = $(wildcard $(DIR)/*.cc) fastsim/xsec.cc
OBJECTS = $(SOURCES:.cc=.o) 
LIB = liblooperBatch.so

DICT = LinkDef_out.o

LIBS = $(LIB)

EXE = main.exe

#
# how to make it
#

$(CORELIB): $(DICT) $(COREOBJECTS)
	$(QUIET) echo "Linking $@"; \
	echo "$(LINKER) -shared -o $@ $(COREOBJECTS) $(DICT) $(LINKERFLAGS)"; \
	$(LINKER) -shared -o $@ $(COREOBJECTS) $(DICT) $(LINKERFLAGS)

$(LIB):	$(DICT) $(OBJECTS) $(COREOBJECTS)
	$(QUIET) echo "Linking $@"; \
	echo "$(LINKER) -shared -o $@ $(OBJECTS) $(COREOBJECTS) $(DICT) $(LINKERFLAGS)"; \
	$(LINKER) -shared -o $@ $(OBJECTS) $(COREOBJECTS) $(DICT) $(LINKERFLAGS)

LinkDef_out.cxx: LinkDef.h
	$(QUIET) echo "Making CINT dictionaries"; \
	rootcint -f LinkDef_out.cc -c -p $(DICTINCLUDE)  LinkDef.h; \
	cat LinkDef_out.cc > LinkDef_out.cxx; rm LinkDef_out.cc

%.exe:  $(LIBS)
	$(QUIET) echo "Building $@"; \
	$(CC) -o $@ $(LIBS) $(ROOTLIBS) ${@:.exe=.cc} 

%.o: 	%.cc %.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cc
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o:    %.cxx 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

libs:	$(LIBS)

build:  $(EXE)

# `make test` will say "up to date" because there's a folder called test
# https://stackoverflow.com/questions/3931741/why-does-make-think-the-target-is-up-to-date
.PHONY: test

b: build

opt: all
small: all
debug: all
tar:
	./make_tar.sh
test:
	./test_babymaker.sh

loopclean:
	rm -f \
	*_out.*	 \
	*.o \
	*.*~ \
	$(LIB) \
	main.exe

clean: loopclean
	rm -f \
	$(CORELIB) \
	./$(COREDIR)/*.o \
	./$(COREDIR)/Tools/MT2/*.o \
	./$(COREDIR)/Tools/*.o \
	./$(COREDIR)/Tools/btagsf/*.o \
	./$(COREDIR)/Tools/JetResolution/*.o \

endif
