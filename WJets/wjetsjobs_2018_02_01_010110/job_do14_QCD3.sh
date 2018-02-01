#!/bin/bash
cd /afs/cern.ch/work/a/awisecar/WJetsTreeAnalysis16/CMSSW_5_3_20/src/WJetsTreeAnalysis16/WJets
eval `scramv1 runtime -sh`

root -l -b -q runDYJets.cc\(14,3\) 2>&1

