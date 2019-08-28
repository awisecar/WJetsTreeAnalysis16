#! /usr/bin/env python2
import json
import pickle

#######################################

whichSF = 1 #ID SF
# whichSF = 2 #Iso SF
# whichSF = 3 #Trig SF

#######################################

dir = '/afs/cern.ch/user/a/awisecar/WJetsTreeAnalysis16_lxplus7/CMSSW_7_6_0/src/WJetsTreeAnalysis16/WJets/EfficiencyTables/sourceFiles/'

if (whichSF == 1):
    jsonfile = '2017_IdSF_BCDEF.json'
elif (whichSF == 2):
    jsonfile = '2017_IsoSF_BCDEF.json'
elif (whichSF == 3):
    jsonfile = '2017_SingleLeptonTriggerSF_BCDEF.pkl'

f = open(dir+jsonfile, 'r')
print '\n---> Opening: '+jsonfile
if not (whichSF == 3):
    results = json.load(f)
else: 
    results = pickle.load(f)

# print(str(results.keys())+'\n')
# print(str(results["IsoMu27_PtEtaBins"].keys())+'\n')

if (whichSF == 1):
    string1 = "NUM_TightID_DEN_genTracks"
    string2 = "abseta_pt"
elif (whichSF == 2):
    string1 = "NUM_TightRelIso_DEN_TightIDandIPCut"
    string2 = "abseta_pt"
elif (whichSF == 3):
    string1 = "IsoMu27_PtEtaBins"
    string2 = "abseta_pt_DATA"

print '\n---> Using: '+string1+', '+string2+'\n'

for etaKey, values in sorted(results[string1][string2].iteritems()):
    for ptKey, result in sorted(values.iteritems()):
        # print "|eta| bin: %s  pT bin: %s\t\tdata/MC SF: %f +/- %f" % (etaKey, ptKey, result["value"], result["error"])
        # print "%s\t%s\t%f\t%f\t%f" % (etaKey, ptKey, result["value"], result["error"], result["error"])
        print "%s\t%s\t\t%.12f\t%.12f\t%.12f" % (etaKey, ptKey, result["value"], result["error"], result["error"])




