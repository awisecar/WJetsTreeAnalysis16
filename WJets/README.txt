# WJetsTreeAnalysis16

# -----------------------------------------------------------------------------------------

purpose --
- this code runs over the events (of both data and MC type) contained in the "Bonzai" trees produced using the ntuple-izer code, performing the object/event selections and filling events that pass these criteria into various histograms for a variety of W+jets observables
- this is done for the data events, and then for the signal MC and BG MC events
- plots are made to assess the agreement between data and predicted (signal+BG) events at reconstructed-level 

set up --
mkdir WJetsTreeAnalysis16_lxplus7
cd WJetsTreeAnalysis16_lxplus7
cmsrel CMSSW_7_6_0
cd CMSSW_7_6_0/src
<pull repo>



running the code --
(running the code in any case first involves the compiling of the code in SourceFiles)

runEvtSelection.cc     >>> the main script that executes the event selection loop for a given set of switches, which are:
  - doWhat - the specific sample to run (the data sample, or the signal MC or one of the BG MC samples)
  - doQCD - selects the W+jets signal region (0) or one of the 3 control regions (1,2,3) used in the data-driven derivation of the QCD BG
  - doSysRunning - whether to run the nominal event selection cuts (0) or the systematic uncertainty variations (all other values), where some of these systematic variations may only run on a subset of the samples used in the analysis
  - doBJets - switch used to include a veto on events that contain >=1 b-tagged jets (-1), or to require events to have >=2 b-tagged jets (2, used for running ttBar control region derive ttBar SFs), or to ignore the b-tagged jet content of the event (0, this is currently the value used for running the analysis with the nominal selections)
  - year - the year of data-taking of the data/MC samples (2016, 2017, 2018 for each of the three years of Run 2)

wjets_jobsub_Condor.py >>> this script submits separate HTcondor jobs for each of the physics processes/QCD region variations/systematic variations/etc needed to run the analysis

--> if you want to run the code locally, comment out the top line in the runEventSelection.cc script and enable the switches doWhat, doQCD, doSysRunning, year --> then run 'root -b -q -l runEventSelection.cc'
--> if you want to run the condor job submission, comment out the mentioned switches and uncomment the top line of the script, then run 'python wjets_jobsub_Condor.py'


in order to add and fill new histograms, need to edit:
SourceFiles/HistoSet.h     >>> declare histogram objects
SourceFiles/HistoSet.cc    >>> give histograms name, title, binning
SourceFiles/ZJetsAndDPS.cc >>> main event loop that makes object/event selection cuts, fills histograms


*** NOTE ***
final reco-level histogram files and plots are already made and are here:
/eos/cms/store/group/phys_smp/AnalysisFramework/Baobab/awisecar/wjetsRun2_histoFiles_forRecoPlots



merging of W+jets histos --

--- to find empty Root files in HistoFiles folders ---
find . -type f -size -450k
---------------------------------------------------

first merge signal region/central, no systematics
./wjets_Merge_Data_NoBVeto.sh
	- OR ./wjets_Merge_Data_BVeto.sh 
./wjets_Merge_Wpt_NoBVeto.sh
	- OR ./wjets_Merge_Wpt_BVeto.sh
merge BG's --
(no TT merge for 2016)
	- OR root -b -q -l runMergeTT.cc
root -b -q -l runMergeTTV.cc
root -b -q -l runMergeVV.cc
root -b -q -l runMergeTop.cc
and run QCD --
root -b -q -l runQCD.cc
then plots --
root -b -q -l runPlots.cc


then merge systematics...  
merge data for syst = 2 --
./wjets_Merge_Data_Syst_NoBVeto.sh
merge W+jets for syst = 1, 4, 5, 11 --
./wjets_Merge_Wpt_Syst_NoBVeto.sh
then merge BG's for syst = 1, 3, 5, 11 --
(no TT merge for 2016)
	- OR root -b -q -l runMergeTT.cc
root -b -q -l runMergeTTV.cc
root -b -q -l runMergeVV.cc
root -b -q -l runMergeTop.cc


======


have to rename histogram folders to 
HistoFiles_2016, etc
all before were named
HistoFiles_2016_noBVeto, etc

then run --
root -b -q -l mergeAllYears.cc  
NOTE: ttbar xsec SFs are applied per year here before merging all years to get Run2 combined
- meaning you have to turn the ttBar xsec SFs off in SourceFiles/Plotter.cc before running plots

to make plots, have to rename HistoFiles_Run2 to HistoFiles before plotting

now make plots for merged Run2 --
root -b -q -l runPlots.cc








to run ttBar SF derivation:
...


to run the theoretical uncertainties (scales, PDF, alpha-s) for the W+jets NLO FxFx samples:
...



# -----------------------------------------------------------------------------------------

##################################################### 
### Note: THIS README (below) NEEDS TO BE UPDATED ###
##################################################### 

The instructions for running the W+jets 13TeV analysis code on bonzai ntuples (data and MC) for CMS Run 2 data.

If you are interested in 2015 data/MC, please see Apichart's repository here:
https://github.com/ahortian/WJetsTreeAnalysis/tree/master/WJets

Note: the unfolding of these data is done at a different step, this code produces reco-level plots.
For unfolding, we are currently using a method based on the TUnfold classes available in ROOT:
https://github.com/awisecar/WJetsUnfolding16

-- For now use CMSSW_5_3_20
cmsrel CMSSW_5_3_20  
cd CMSSW_5_3_20/src
cmsenv 

-- Clone /WJetsTreeAnalysis from the github
git clone https://github.com/awisecar/WJetsTreeAnalysis16.git
or
git clone git@github.com:awisecar/WJetsTreeAnalysis16.git

-- Go to the directory 
cd WJetsTreeAnalysis16/WJets

-- Compile RooUnfold
cd RooUnfold
make clean
make
cd -

-- Ready to run the code!!! 

Two different approaches to running:

1) python wjets_jobsub_Condor.py
which will submit jobs to HTCondor for the values of doWhat, doQCD, and doSysRunning specified in the script, producing STDOUT, and the submitted shell scripts in the folder "wjetsCondor_"

2) root -b -q -l wjets_compileCode.cc
   root -b -q -l runDYJets.cc
first compiles the code anew, and then the main command, which will be run in your current shell.

Note: for QCD multijets background estimation we use the "ABCD" data-driven method, so we must run over the one signal and three control regions, doQCD=0, 1...3, respectively.
The doWhat, doQCD, doSysRunning, and doBJets options tell the main code (Sources/ZJetsAndDPS.cc) what you'd like to run

-- Output root files are in the directory HistoFiles. 

-- Merge each of the separate datafiles. (The data batch jobs are split into separate lists -- see runDYJets.cc -- and submitted separately to speed up the process.) This shell script uses the 'hadd' command to combine several rootfiles into one.
./wjets_Merge_Data_BVeto.sh
./wjets_Merge_Wpt_BVeto.sh

-- Merge single top and DYJets root files (change if your total number of top or DY MC samples changes):
root -b -q -l MergeTop_BVeto.cc++
root -b -q -l MergeVV_BVeto.cc++

-- Merge QCD Background:
root -b -q -l runQCD_BVeto.cc++
or if you want to save the output --
root -b -q -l runQCD_BVeto.cc++ > mergeQCD_output.out 2>&1 

-- Plotting: After you have all needed output files, you can plot the histograms by
root -b -q -l runTest_BVeto.cc++

--Congratulations! Output PDF files will be in the directory PNGFile
A .tex file detailing the bin statistics will be available in WJets directory with title beginning with "statsTable_"

----------------------------------

Systematics:
./wjets_Merge_Data_Syst.sh
./wjets_Merge_Wpt_Syst.sh
//Note: for single top, go in and change the syst switches as needed, but it's the same file to merge them
root -b -q -l MergeTop_BVeto.cc++

----------------------------------

TTbar SFs Study:
Run requiring 2 btag jets (doBJets=2) and WJets signal region (doQCD=0)
We neglect QCD BG here.

./wjets_Merge_Data_BJets.sh
./wjets_Merge_Wpt_BJets.sh

root -b -q -l MergeTop_BJets.cc++
root -b -q -l MergeVV_BJets.cc++

//if(doTTScale), then ttbar SFs are applied in getStatistics function to print out the statistics tables
//ttbar SFs are applied in Plotter.cc when plotting the distributions (but original files not altered)
root -b -q -l runTest_BJets.cc++
