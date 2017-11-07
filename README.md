# WJetsTreeAnalysis16



The instructions for running the W+jets 13TeV analysis code on bonzai ntuples (data and MC).

-- For now use CMSSW_5_3_20
cmsrel CMSSW_5_3_20  
cd CMSSW_5_3_20/src
cmsenv 

-- Clone /WJetsTreeAnalysis from the github
git clone https://github.com/awisecar/WJetsTreeAnalysis16.git

-- Go to the directory 
cd WJetsTreeAnalysis16/WJets

-- Compile RooUnfold
cd RooUnfold
make clean
make
cd -

-- Ready to run the code!!! 

-- Run the code under /WJets
root -b runDYJets.cc++

-- For QCD Background go to WJets/runDYJets.cc change int doQCD = 1
root -b -q runDYJets.cc
# Change it 1-3 and run every time. 

-- Output root files are in the directory HistoFiles.

-- Merge single top and DYJets root files:
root -b -q MergeTop_BVeto.cc++
root -b -q MergeVV_BVeto.cc++

-- Merge QCD Background:
root -b -q runQCD_BVeto.cc++

-- Plotting: After you have all needed output files, you can plot the histograms by
root -b -q runTest_BVeto.cc++

--Congratulations! Output PDF files will be in the directory PNGFiles
