# Makefile

CC = g++

NTUINC    = -I/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/tHFCNC/NtupleProducer/
NTULIB    = -L/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/tHFCNC/NtupleProducer/\
              -lNtuple
	      
JECLIB      = -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_12/lib/slc6_amd64_gcc530/ -lCondFormatsJetMETObjects -lJetMETCorrectionsModules -lRecoEgammaEgammaTools

SOVERSION               = $(shell echo $(VERSION) | \
                               sed 's/\(.*\)\.\(.*\)\.\(.*\)/\1:\2:0/')
ROOTCONF                = root-config
ROOTCINT                = $(shell $(ROOTCONF) --prefix)/bin/rootcint
ROOTLIBS                = $(shell $(ROOTCONF) --libs)
ROOTGLIBS               = $(shell $(ROOTCONF) --glibs)
ROOTCFLAGS              = $(shell $(ROOTCONF) --cflags)
INCLUDES                = $(shell $(ROOTCONF) --incdir) $(ROOTCFLAGS)
CINTINCLUDES            = $(shell $(ROOTCONF) --incdir)
LIBS                    = $(ROOTLIBS)

LDLIBS  = $(ROOTLIBS) $(NTULIB) $(JECLIB)

ADDINC = $(NTUINC)

SRC = gen.cxx

INC = 

OBJ = gen.o

EXE = gen

all:  	     		$(EXE)
#			@rm -f *.o

$(EXE):     	 	$(OBJ) $@

$(OBJ):     		$(SRC) $(INC)
		        $(CC) -c $(ROOTCFLAGS) $(ADDINC) $(SRC)

clean: 		
			@rm -f $(OBJ) $(EXE)
