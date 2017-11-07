# LSBATCH: User input  

USER_DIR="/afs/cern.ch/user/a/awisecar"
CMSSW_PROJECT_SRC="WJetsTreeAnalysis16/CMSSW_5_3_20/src"
SELECTOR_DIR="WJetsTreeAnalysis16/WJets"

cd $USER_DIR/$CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
cd $USER_DIR/$CMSSW_PROJECT_SRC/$SELECTOR_DIR
root -b -q runDYJets.cc