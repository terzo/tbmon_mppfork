CXX          = g++
CXXFLAGS     = -g
LD           = g++
LDFLAGS      = -g -L.

ROOTCONFIG   := $(shell which root-config)
ROOTCINT     := $(shell which rootcint)

ROOTCXXFLAGS := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
LIBS         := $(USERLIBS) $(ROOTLIBS) $(LIBS)

INCLUDEDIRS  := -
CXXFLAGS     := $(ROOTCXXFLAGS) $(CXXFLAGS) -fPIC 
LDFLAGS      := $(ROOTLDFLAGS) $(LDFLAGS) 


all: totcal

totcal:	totcal.cpp rootlogonCompile.h
	$(LD) $(LDFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS) 

clean:
	rm *.o totcal
