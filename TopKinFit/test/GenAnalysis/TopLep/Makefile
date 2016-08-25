# Makefile

CC = g++

NTUINC    = -I/home-pbs/kskovpen/ttH/CMSSW_7_6_3_patch1/src/IPHCNtuple/NtupleProducer/
NTULIB    = -L/home-pbs/kskovpen/ttH/CMSSW_7_6_3_patch1/src/IPHCNtuple/NtupleProducer/\
              -lNtuple

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

LDLIBS  = $(ROOTLIBS) $(NTULIB)

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
