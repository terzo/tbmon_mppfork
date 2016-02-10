/** @mainpage 

User/developer manual for the ATLAS 3D pixel analysis framework.

@author Havard Gjersdal
@author Kyrre Sjøbæk (cmdLineExtras)


This document tries to explain how to use the ATLAS 3D pixel analysis framework, as well as how to expand it. 


@section introduction Introduction
Tbmon is a small analysis framework for the ATLAS 3D pixel test beams built to study the output root
files from the track reconstruction.
For a given set of runs, it produces a set of plots for all specified devices under test.
The framework is intended to be flexible enough to run the same analysis on data from different beam
telescopes, and modular enough for people to work on different parts of the analysis without
stepping on each others toes

@section quick Quick Start

Get some data:
@code
$ mkdir -p ~/tb2009/data
$ cd ~/tb2009/data
# setup ssh tunnel, replace <USERNAME> with your cern login
$ ssh -f -L 8000:mpnoatlas01.cern.ch:8000 <USERNAME>@@lxplus.cern.ch -N
# download all root files on server for may2009(requires wget)
$ wget -r -nd -A "*.root" http://localhost:8000/~atlas3dsi/tb2009/tbtrack/
@endcode

Get the code:
@code
$ cd ~/tb2009
# replace <USERNAME> with your cern login
$ svn co svn+ssh://<USERNAME>@@svn.cern.ch/reps/atlas3dpix/tbonline
@endcode

In order to get an analysis running, it needs to know where to find data and where to place it's output. 
These configurations are expected to be found in siteconfig.h. This file is not under version
control, so copy siteconfig.h.example  to siteconfig.h and set up the paths correctly.
@code
trySet(config.outPath, (char*) "~/tb2009/output/");
trySet(config.plotExtension, (char*) ".pdf");
trySet(config.dataPath, (char*) "~/tb2009/tbtrack/");
@endcode
The paths need to exist. If the extension is not set, or set to something root does not reckognize,
it will default to eps. Stuff set in the siteConfig function will override the other functions, so
if you have data from may on place, and data from october somewhere else, do not set datapath in
siteConfig, but use the config spesific functions.

Compile:
@code
make
@endcode

Running ./tbmon should produce a full and updated list of available command line arguments. Something like
@code
./tbmon -a residuals -i 160 -l runlists/may2009-bon-a15
@endcode
will run  the analysis named "residuals" on DUT with identifyer 160 over all runs in the runlist
runlists/may2009-bon-a15.

@section compile Compiling/running


In order to build the program, ROOT must be installed on the system, rootcint and root-config must be
found in $PATH.

To build, simply:
@code
$ make
@endcode

Some older versions of ROOT can not open the tbtrack#.root files without
a proper description of the Track class. This can be fixed by compiling
the Track class and loading it into ROOT before opening the file.
@code
$ make libTrack.so
$ make libModule.so
$ root -l
root [0] .L libTrack.so
root [1] .L libModule.so
root [2] TFile* f = new TFile("tbtrackXXXX.root")
root [3] TBrowser tb
@endcode

To rebuild this documentation:
@code
make docs
@endcode

@section running Running
Run ./tbmon for a full and updated list of command line arguments. When this was written:
@code
Needs ONE of the folowing arguments to determine what files to loop over:
 -l <runList>                    Loop over all runs in <runList>
 -r <firstrun> <lastrun>         Loops for all runf from <firstrun> to <lastrun>
 -s <run>                        Loop over a single run

Optional arguments:
 -a <analysisname>               only run analysis named <analysisname>. If not supplied it runs over all available analyses.
 -i <iden>                       only run over iden <iden>. If not suppliedit loops over all available idens.
 -c <configname>                 Configures modules for configuration <configname>. If not supplied, by default defaults to may2009
 -v <verbosity-level>            Should be one of (in order of increasing verbosity) error, info, debug, debug2 or debug3. Default is info.
 -e <extension>                  Plot extension. Should be a file format reckognized by ROOT.
 -d <data-path>                  Path to input tbtrack????.root files.
 -o <out-path>                   Path to output directory.
 -b <base-name>                  Base name for all output. Defaults to name based on run list or range of runs.
 -f                              Redirect output(stdout) to log file named after the base name.
 -P:<key> <value>                Pass an extra parameter, which can be used by analyses or builders

Simple Example:
 tbmon -l runlists/may2009-bon-a15
Applies all available analyses to all the runs in runlist for all idens

Complete Example:
 tbmon -l runlists/may2009-bon-a15 -i 160 164 -a sumtot residuals -c may2009 -b test1 -e png -f -o out -P:param1 1.3 -P:param2 urk
Applies analysis named "sumtot" and analysis named "residuals" to all runs in runlist for idens 160 and 164. Png plots are made and saved to relative path out. A log file(out/test1.log) is made, no output to screen. Extra-parameters "param1" and "param2" is stored internally with values "1.3" and "urk", and may be accessed from any builder and/or analysis.
@endcode

@subsection cmdLineExtras ``-P:<key> <value>''arguments

These special arguments are usable for sending special options to eventbuilders or analyses.
For example, to set the chi2 cut to 10.0 through the chi2builder class, use the option
@code
-P:G_chi2cut 10.0
@endcode
All such parameters currently defined is printed togeter with their current value
and a short help-text before starting to process runs.
@code
[ TbConfig ]; Parameter arguments pushed onto cmdLineExtras:
        KEY                           VALUE               DESCRIPTION
        B_simDutRunner_model          Full3D_HP           Name of model wanted. Leave to default ("passthrough") to use PllHits from .bdt
        G_chi2cut                     10                  Mark tracks with chi2 > this value as bad
        SD_Full3D_HP_R_bias           5                   Radius of bias electrodes    [um]
        SD_Full3D_HP_R_readout        5                   Radius of readout electrodes [um]
@endcode

The implementation of this, and how to use it in your analysis/eventbuilder is described
in the configuration section.


@section configuration Configuration

Configuration is done by manipulating an instance of the TbConfig class. This can be done by passing
in command line arguments(see Compiling/running), or by editing the driver.cc/.h/siteconfig.h
files. The command line arguments should allways override everything else.

@section site Site configuration
In order to not litter the repository with site spesific configurations, this should be done in the
file siteconfig.h, which is not under version control. A good place to start setting this up is by
copying siteconfig.h.example, which is in version control. In addition to setting up input/output
paths, the defaults of any part of TbConfig can be set here, like the default verbosity level, the
default idens or analyses to loop over.

@section calibration Calibrations
There are two intended ways on storing calibrations in the framework. Configurations that are
specific for a DUT is stored in a DUT object. This will include stuff lie mask, eta correction and
tot to charge conversion calibrations. Global calibrations, like track angle cut, TDC, and chi2 cut
calibrations should be stored in objects inheriting from the EventBuilder class. Both these classes
have an initRun function that makes it possible to extract the right calibration for a spesific run,
if this is needed. The DUT's and EventBuilders are in turn stored in the TbConfig object.

These calibrations should always be read in from text/root files, and never hard coded.

Since many of these calibrations will vary between the different test beam configurations, the plan
is to read in the calibrations in one function per test beam configuration in driver.cc.

If you update any of these calibrations or cut values, be sure to announce it on the mailing list
when you commit.

@section extras cmdLineExtras

As discussed in the Compiling/running section, there exists a system for passing ``special''
flags and arguments to builders or analyses. This may be where to set cuts, histogram binning/ranges,
name of wanted simdut model etc.

The cmdLineExtras mechanism is implemented through two string-indexed hashmaps of strings in
TbConfig: cmdLineExtras and cmdLineExtras_help. To use this in your own analysis etc., you have to
add code setting defaults and parsing it in your init() method. As an example,
see chi2builder:
@code
cutVal = config.cmdLineExtras_argGetter("G_chi2cut", chi2default,"Mark tracks with chi2 > this value as bad");
@endcode
This uses the method
@code
double TbConfig::cmdLineExtras_argGetter(string key, double default, string description)
@endcode
to create a new key (a new parameter) and set its default value and description.
It then returns either the default value, or a value previously stored (for example,
from the command-line, or from the siteconfig via TbConfig::cmdLineExtras_trySet() ).
As cmdLineExtras_argGetter is overloaded, there is another such method for strings:
@code
string TbConfig::cmdLineExtras_argGetter(string key, string defaultValue, string describeText)
@endcode
You are more than welcome to contribute methods for other types, ranges (returning vectors), etc.

@subsection naming Naming convention

As cmdLineExtras is a ``namespace'' shared by all parts of the code, it is important to not mess it up.
Thus I propose a naming scheme for the keys:
@itemize
@item G_<name>: Global things, used by several parts of the framework
@item B_<buildername>_<varname>:  Arguments to the builder <buildername>
@item A_<analysisname>_<varname>: Arguments to the analysis <analysisname>
@item SD_<sdname>_<varname>:      Argument to the simDut <sdname>
@end itemize


@section core The core of the framework
@section tbconfig TbConfig
The TbConfig contains all information about the running job. In addition to keeping track of all
running analysis classes and all calibrations/configurations it contains methods for writing
histograms to image files or a root file.
@section event The Event
The Event object contains all information about the current event. It is created for the DUT that
the current TbAnalysis object is studying, and contains a pointer to this DUT object.
@itemize
@item event.track: A pointer to the current track object. This object contains all information about the current event(Note: Only available for BAT data!).
@item event.hits: This is a vector of all PllHits* belonging to the current DUT
@item event.clusters: This is a vector of clusters belonging to the current DUT. Clusters are vectors of PllHits*, so it
is a vector of vectors.
@item event.trackX, event.trackY, event.dxdz, event.dydz: The fitted track parameters extrapolated to the current DUT.
@end itemize
@subsection flags Flags
The event object also contains a set of flags, which are intended to be used for cuts in the
analysis objects.
@itemize
@item fTrack: Should equal event::kGood if the track has passed all cuts.
@item fTrackChi2: Equals event::kGood if the track passed chi2 cuts.
@item fTrackAngle: Equals event::kGood if the track passed angle cuts.
@item fTrackCentalRegion: Equals event::kGood if track passes through central region of DUT. This
does not affect the fTrack flag.
@item fClusters: Equals event::kGood if tracks are built successfully. Cluster quality checks are
not yet implemented.
@item fEtaCorrections: Equals event::kGood if eta corrections are available.
@end itemize
@subsection constants Constants
The event namespace also contains a set of constants to be used in the analysis
@itemize
@item  event::pitchX = 400.0;
@item  event::pitchY =  50.0;
@item  event::ncols = 18;
@item  event::nrows = 160;
@item  event::maxCol = 17;
@item  event::maxRow = 159;
@item  event::skipRows = 16;
@item  event::skipCols = 2;
@item  event::matchX = 400;
@item  event::matchY = 150;
@end itemize
@section eventbuilder The EventBuilders
The events are built by classes inheriting from the EventBuilder class. Examples are cluster
finders, that fills the event.clusters vector, or track quality checkers that set the different
track flags.
@section cluster Cluster
Clusters are not stored in a special class, and are just vectors pf PllHits*. Methods for dealing
with clusters are available in the cluster namespace. See core/include/clusters.h for details.
@itemize
@item int    getSumTot( const event::cluster_t &cluster ): Returns the cluster ToT 
@item int    getMaxTotCluster(const vector< event::cluster_t > &clusters): Returns the vector index
of the maxTot cluster
@item int    getMatched(const vector< event::cluster_t > &clusters, const Event &event): Returns the
  vector index of the cluster matching the track hit position. -1 if no match.
@item int    getMaxTotCell(const event::cluster_t &clusters): Returns the vector index of the max tot cell
in a cluster.
@item double getUnWeightedX(const event::cluster_t &cluster): Geometrical mean of the cluster
@item double getUnWeightedY(const event::cluster_t &cluster): 
@item double getChargeWeightedX(const event::cluster_t &cluster): Charge weighted mean of the cluster
@item double getChargeWeightedY(const event::cluster_t &cluster):
@item double getEtaCorrectedX(const event::cluster_t &cluster, DUT* dut): Corrected charge weighted mean
of the cluster. This is the most accureat position estimator.
@item double getEtaCorrectedY(const event::cluster_t &cluster, DUT* dut):
@end itemize
@section pixelhit pixelhit
The PllHit class does not contain any methods, so a few have been implemented in the pixelhit namespace

@itemize
@item bool inCentralRegion(PllHit* hit): Is the hit in the central region of the DUT?
@item bool isMatch(PllHit* hit, const Event &event): Does the hit match the extrapolated track?
@item int getCol(double trackX): What column number does trackX pass through? Returns -1 if out of
bounds. Does not take active edges into account.
@item int getRow(double trackY): What row number does trackY pass through? Returns -1 if out of
bounds. Does not take active edges into account. 
@end itemize

@section logging Logging
A primitive logging system has been implemented using only cout. Passing the -f flag to tbmon will
redirect the standard output to a file named after the base name of the analysis. If you want text
to both logs and  screen something like
@code
./tbmon <arguments> |tee logfile.log 
@endcode
should do the trick.

The TbConfig object contains a verbosity level, set to one of kERROR, kINFO, kDEBUG, kDEBUG2,
kDEBUG3. 

A pre-compiler macro is available for TbAnalysis classes, this should be used instead of
calling cout directly:
@code
TBALOG(kDEBUG3) << "Found matching cluster!" << endl;
@endcode
will expand to
@code
if( kDEBUG3 <= config.logLevel) cout << "[ " << this->name << " ]: " << "Found matching cluster!" << endl;
@endcode
and print something like 
@code
[ efficiency-160 ]: Fount matching cluster!
@endcode
as long as the verbosity level is higher or equal to kDBUG3.

Please use kDEBUG2 or kDEBUG3 for stuff that is printed every event.

@section add Adding analysis code



The analysis is performed by objects inheriting from the TbAnalysis class.

@section guidelines Guidelines
Several people might be interested in using/improving on your code. It is therefor important to try
and keep it readable by others.

@itemize
@item The analysis class should implement a(one) specific analysis. If the same class plots hit
efficiency and charge sharing, you are most definitely doing it wrong.
@item Try to keep the number of produced plots to a reasonable minimum. Plots for debugging should
not be printed when debugging is done.
@item Analysis classes should be very quiet at info and error verbosity levels. Use the TBALOG macro. 
@item If you end up copying and pasting methods between different analysis classes, this means it
probably belongs somewhere in the core framework.
@item Do not commit code that does not run, but other than that commit often. People might be
interested in what you are doing.
@item After test beams and before big presentations  the code will be actively developed. Update often.
@item If you are having problems, look for help in other analysis classes. Use the mailing list:
atlas-project-3D-software@@cern.ch 
@end itemize

@section virtual Virtual functions

@itemize
@item void init(TbConfig &config): 
This function is called before the loop over root files. Histograms, member variables should be initialized here.
Purely virtual, must be implemented.
@item void initRun(const TbConfig &config):
This function is called before each root file is processed. In most cases this function does not need to be implemented.
@item void event(const TbConfig &config, const Event &event):
This function is called for each entry in the root files. This is where the main part of the
analysis is done. Purely virtual, must be implemented.
@item void finalizeRun(const TbConfig& config):
This function is called after each root file has been processed. In most cases this function does not need to be implemented.
@item void finalize(const TbConfig &config):
This function is called when all the root files are processed. Histograms should be saved to files, efficiencies or similar
calculated and printed. Purely virtual, must be implemented.
@end itemize 

@section driver driver.cc and Makefile
To include the analysis in the framework, you need to do a few things:
@itemize
@item Source files should have extensions .cc and .h to keep Makefile simple
@item Add <filename>.o to ANALYSIS list in Makefile
@item Add header to driver.cc
@item Add class to framework in allAnalyses function in driver.cc.
@item The name you give it in allAnalyses is the name you call it by from the command line.
@end itemize

@section obtain_calibration Obtaining calibrations
Most of the calibrations can be obtained by running an analysis job that produces a text file, that
can be read back in to the driver.cc file.
@itemize
@item checkdutsync: Checks for windows where the sync between the DUT and the tracks are
uncertain. Should be run over all runs for configuration. Only needs one DUT enabeled.
@code
./tbmon -c may2009 -a checkdutsync -r 500 1500 -i 160
@endcode
@item angledist: Checks the mean and rmsd of the fitted track angles. Used for angle cuts. Should
run over all runs in a configuration, only one DUT needs to be enabeled.
@code
./tbmon -c may2009 -a angledist -r 500 1500 -i 160
@endcode
@item getetacorr: Gets the eta corrections for a run. This must be run over all DUTs, one file is
needed for each change in the setup that will affect the eta correction(changed tilt/ magnetic field).
@code
./tbmon -c may2009 -a getetacorr -l runlists/may2009-bon-a15
@endcode
@end itemize

Calibrations for the tot->q conversion can not be obtained by running an analysis job. Calibrations
for masks have not yet been implemented as an analysis job.

@section about About this document
This document is created with texinfo. Feel free to extend it by editing the tbmon.texi file. To
rebuild the info and html files: 
 
@code
make docs
@endcode

If there are any questions or comments, send a mail to the mailing list:

atlas-project-3D-software@@cern.ch 



*/
