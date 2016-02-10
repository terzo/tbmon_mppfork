CXX          = g++
CXXFLAGS     = 
LD           = g++
LDFLAGS      = -L.

ROOTCONFIG   := $(shell which root-config)
ROOTCINT     := $(shell which rootcint)

ROOTCXXFLAGS := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
LIBS         := $(USERLIBS) $(ROOTLIBS) $(LIBS)

INCLUDEDIRS  := -
CXXFLAGS     := $(ROOTCXXFLAGS) $(CXXFLAGS) -fPIC 
LDFLAGS      := $(ROOTLDFLAGS) $(LDFLAGS) 

CORE          = checktrack.o clusters.o dut.o totcalib_fei3_turbodaq.o totcalib_fei3_usbpix_converted.o etacorrections.o looper.o totcalib_fei3_turbodaq_converted_from_tot_file.o simThreeVector.o tbconfig.o tbutils.o totcalib.o Track.o TrackDict.o

EVENTBUILDERS =  anglecuts.o battrack.o calcangles.o checkregion.o  checkcentralregion.o chi2builder.o clusterdumper.o clusterfinder.o clustermasker.o dutsync.o batetacutter.o eubuildtrack.o maskandlvl1.o maskreader.o lvl1cuts.o pixelmasker.o totcalibreader.o translator.o translatorRunningXY.o simBaseBuilder.o simDutRunner.o simPixelEdepBuilder.o simTruthBuilder.o

ANALYSIS      = batangledist.o batunbiased.o beamprofile.o blank.o checkalign.o checkalignRunningXY.o batcheckdutsync.o checktrack.o clusterchecker.o clusterpixtot.o clustersvsrun.o correlations.o edgeefficiency.o edgeefficiencyshift.o efficiency.o efficiency2.o efficiencysimp.o efficiencyvsrun.o etawidth.o getetacorr.o hotpixelfinder.o lvl1cut.o maxcellresiduals.o qEfficiency.o qshare1D.o qshare2D.o readout.o residuals.o sumtot.o simDutEdep.o simResiduals.o botho.o

SIMDUT        = pixel_simple.o Full3D_HP.o Full3D_Vadim.o

STYLE = AtlasStyle.o AtlasLabels.o

OBJS = $(addprefix temp/,$(CORE) $(ANALYSIS) $(EVENTBUILDERS) $(SIMDUT) $(STYLE))

vpath %.cc core/src analysis/src eventbuilders/src simdut/src style/src
vpath %.h core/include analysis/include eventbuilders/include simdut/include style/include

DEBUG ?= 1
ifeq (DEBUG, 1)
    CXXFLAGS +=-DDEBUG -g
    LDFLAGS += -g
else
    CXXFLAGS +=-DNDEBUG -O2
endif

all: tbmon setup

#Force everything to recompile if event.h updates.
#-O3 -msse2 -ftree-vectorize -ftree-vectorizer-verbose=5
temp/%.o : %.cc %.h event.h
	$(CC) -c -Icore/include/ -Ianalysis/include -Ieventbuilders/include -Isimdut/include -Istyle/include $(CXXFLAGS) $< -o $@

libTrack.so: temp/Track.o temp/TrackDict.o
	$(LD) $(LDFLAGS) -Icore/include/ -shared -o $@ $^

TrackDict.cc: Track.h TrackLinkDef.h
	$(ROOTCINT) -v4 -f $@ -c $^

tbmon: driver.cc driver.h siteconfig.h $(OBJS)
	mkdir -p temp
	$(LD) $(LDFLAGS) -Icore/include/ -Ianalysis/include -Ieventbuilders/include -Isimdut/include -Istyle/include $(CXXFLAGS) -o $@ $^ $(LIBS)
	
setup: setupdriver.cc driver.h siteconfig.h $(OBJS)
	mkdir -p temp
	$(LD) $(LDFLAGS) -Icore/include/ -Ianalysis/include -Ieventbuilders/include -Isimdut/include -Istyle/include $(CXXFLAGS) -o $@ $^ $(LIBS)

docs: doc/Doxyfile
	cd doc && doxygen

#sql23: sql23.cc temp/Track.o temp/TrackDict.o temp/module.o # Doesn't work...
#	$(LD) $(LDFLAGS) -Icore/include/ -o $@ $^ $(LIBS)


clean:
	rm temp/* TrackDict.* tbmon
	
clean-doc:
	rm -rf doc/html/* doc/latex/* doc/man/*

fitter: APIXFitter.cc APIXFitter.h
	$(LD) $(LDFLAGS) $(LIBS) -lgsl -lgslcblas -Icore/include/ -Ianalysis/include -Ieventbuilders/include $(CXXFLAGS) -o $@ $^ 
