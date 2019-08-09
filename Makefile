PROGRAM       = JEWELChargedDijet

version       = CDijet
CXX           = g++
CXXFLAGS      = -O -Wall -g -Wno-deprecated -D$(version)
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
#############################################
# -bind_at_load helps to remove linker error
############################################
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS  = $(shell root-config --evelibs)
CXXFLAGS += $(shell $(FASTJET)/bin/fastjet-config --cxxflags )
LDFLAGS += $(shell $(FASTJET)/bin/fastjet-config --libs --plugins ) 
LDFLAGS += -L$(HEPMC)/lib -lHepMC 
LDFLAGS += -Wl,-rpath -Wl,$(HEPMC)/lib
LDFLAGS += -L$(HEPPDT)/lib -lHepPDT -lHepPID
LDFLAGS += -Wl,-rpath -Wl,$(HEPPDT)/lib
INCS    += -I$(HEPMC)/include
INCS    += -I$(HEPPDT)/include
LDFLAGS += -L/home/alidock/.sw/slc7_x86-64/cgal/v4.6.3-48/lib/ -lCGAL 
LDFLAGS += -L/home/alidock/.sw/slc7_x86-64/GMP/v6.0.0-45/lib/ -lgmp
CXXFLAGS  += $(INCS)
LDFLAGS += $L -ldl

HDRSDICT = src/AliJCDijetHistos.h \
           src/AliJCDijetAna.h \
		   src/AliJHistogramInterface.h \
		   src/AliJHistManager.h \
		   src/AliJBaseTrack.h \
		   src/AliJBaseCard.h \
		   src/AliJBaseEventHeader.h \
		   src/JTreeDataManager.h \
		   src/iaaAnalysis/AliJIaaAna.h \
		   src/iaaAnalysis/AliJIaaHistograms.h \
		   src/iaaAnalysis/AliJIaaCorrelations.h \
		   src/AliJCorrelationInterface.h \
		   src/AliJPiZero.h \
		   src/AliJPhoton.h \
		   src/AliJMCTrack.h \
		   src/AliJTrackCut.h \
		   src/AliJEventHeader.h \
		   src/AliJRunHeader.h \
		   src/AliJDataManager.h \
		   src/AliJCard.h \
		   src/AliJEventPool.h \
		   src/AliJEfficiency.h \
		   src/AliJRunTable.h \
		   src/AliJTrack.h 

HDRS	+= $(HDRSDICT)  nanoDict.h


SRCS = $(HDRS:.h=.cxx)
OBJS = $(HDRS:.h=.o)

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS) src/AliJConst.h $(PROGRAM).C
		@echo "Linking $(PROGRAM) ..."
		$(CXX)  -lPhysics -L$(PWD) $(PROGRAM).C $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $(PROGRAM) 
		@echo "finally done"

%.cxx:

%: %.cxx
#  commands to execute (built-in):
	$(LINK.cc) $^ $(CXXFLAGS) $(LOADLIBES) $(LDLIBS) -o $@

%.o: %.cxx %.h
#  commands to execute (built-in):
	$(COMPILE.cc) $(OUTPUT_OPTION) $<


clean:
		rm -f $(OBJS) core *Dict* $(PROGRAM).o *.d $(PROGRAM) $(PROGRAM).sl

cl:  clean $(PROGRAM)

nanoDict.cc: $(HDRSDICT)
		@echo "Generating dictionary ..."
		@rm -f nanoDict.cc nanoDict.hh nanoDict.h
		@rootcint nanoDict.cc -c -D$(version) $(HDRSDICT)
