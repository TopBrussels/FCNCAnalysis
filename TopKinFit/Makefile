ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lMinuit
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = gcc -fPIC
LD            = gcc
LDFLAGS       = -g
SOFLAGS       = -shared

CXXFLAGS       = $(ROOTCFLAGS)
INCLUDE_FLAGS  = 
LDLIBS         = $(ROOTLIBS)
GLIBS          = $(ROOTGLIBS)

EXE           = kinfit

INC 	      = kinfit.h TopTopLepLep.h TopLep.h TopTopLepHad.h \
TopTopLepHbb.h TopHLepbb.h

SRC	      = kinfit.cxx TopTopLepLep.cxx TopLep.cxx \
TopTopLepHad.cxx TopTopLepHbb.cxx TopHLepbb.cxx

OBJS          = kinfit.o TopTopLepLep.o TopLep.o TopTopLepHad.o \
TopTopLepHbb.o TopHLepbb.o

LIB           = libKinFit.so

all: 	      $(LIB)

$(LIB):	      $(INC) $(SRC)
	      @echo "####### Generating dictionary"
	      @rootcint -f kinfitDict.cxx -c -p $(CXXFLAGS) \
	      $(INCLUDE_FLAGS) -I. $(INC) LinkDef.h

	      @echo "####### Building library $(LIB)"
	      @$(CXX) $(SOFLAGS) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDE_FLAGS) -I. $(SRC) \
	      kinfitDict.cxx -o $(LIB) $(ROOTLIBS)
	      
	      @echo  "####### Removing generated dictionary"
	      @rm -f kinfitDict.cxx kinfitDict.h
	      @rm -f *.o

clean:
	      @rm -f $(OBJS) $(EXE) kinfitDict.cxx kinfitDict.h $(LIB)
