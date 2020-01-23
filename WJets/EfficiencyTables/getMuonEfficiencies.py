#! /usr/bin/env python2
import json
import pickle

#######################################

# year = 2016
# year = 2017
year = 2018

# whichSF = 1 #ID/Tracks SF
whichSF = 2 #Iso/ID SF
# whichSF = 3 #Trig/Iso SF

#######################################

dir = '/afs/cern.ch/user/a/awisecar/WJetsTreeAnalysis16_lxplus7/CMSSW_7_6_0/src/WJetsTreeAnalysis16/WJets/EfficiencyTables/muonEffSFs/'

if year == 2016:    ### currently divided into eras
    if (whichSF == 1):
        jsonfile = '2016Legacy_IdSF_BCDEF.json'
        # jsonfile = '2016Legacy_IdSF_GH.json'
    elif (whichSF == 2):
        jsonfile = '2016Legacy_IsoSF_BCDEF.json'
        # jsonfile = '2016Legacy_IsoSF_GH.json'
    elif (whichSF == 3):         
        jsonfile = ''
        print(" >>>>> No Trig/Iso SFs for 2016 legacy rereco yet! <<<<<")
elif year == 2017:
    if (whichSF == 1):
        jsonfile = '2017_IdSF_BCDEF.json'
    elif (whichSF == 2):
        jsonfile = '2017_IsoSF_BCDEF.json'
    elif (whichSF == 3):
        jsonfile = '2017_SingleLeptonTriggerSF_BCDEF.pkl'
elif year == 2018: 
    if (whichSF == 1):
        jsonfile = '2018_IdSF_ABCD.json'
    elif (whichSF == 2):
        jsonfile = '2018_IsoSF_ABCD.json'
    elif (whichSF == 3):
        jsonfile = ''

f = open(dir+jsonfile, 'r')
print '\nOpening: '+jsonfile
if not (whichSF == 3):
    results = json.load(f)
else: 
    results = pickle.load(f)

# print(str(results.keys())+'\n')
# print(str(results["IsoMu27_PtEtaBins"].keys())+'\n')

if year == 2016:
    if (whichSF == 1):
        string1 = "NUM_TightID_DEN_genTracks"
        string2 = "pt_eta"
    elif (whichSF == 2):
        string1 = "NUM_TightRelIso_DEN_TightIDandIPCut"
        string2 = "pt_eta"
    elif (whichSF == 3):
        string1 = ""
        string2 = ""
elif year == 2017:
    if (whichSF == 1):
        string1 = "NUM_TightID_DEN_genTracks"
        string2 = "abseta_pt"
    elif (whichSF == 2):
        string1 = "NUM_TightRelIso_DEN_TightIDandIPCut"
        string2 = "abseta_pt"
    elif (whichSF == 3):
        string1 = "IsoMu27_PtEtaBins"
        string2 = "abseta_pt_DATA"
elif year == 2018:  
    if (whichSF == 1):
        string1 = "NUM_TightID_DEN_TrackerMuons"
        string2 = "abseta_pt"
    elif (whichSF == 2):
        string1 = "NUM_TightRelIso_DEN_TightIDandIPCut"
        string2 = "abseta_pt"
    elif (whichSF == 3):
        string1 = ""
        string2 = ""

print 'Using: '+string1+', '+string2+'\n'

if year == 2016:
    for ptKey, values in sorted(results[string1][string2].iteritems()):
        for etaKey, result in sorted(values.iteritems()):
            # print "%s\t%s\t\t%.12f\t%.12f\t%.12f" % (ptKey, etaKey, result["value"], result["error"], result["error"])
            print "%s\t%s\t\t%.12f\t%.12f\t%.12f" % (etaKey, ptKey, result["value"], result["error"], result["error"])
elif year == 2017:
    for etaKey, values in sorted(results[string1][string2].iteritems()):
        for ptKey, result in sorted(values.iteritems()):
            print "%s\t%s\t\t%.12f\t%.12f\t%.12f" % (etaKey, ptKey, result["value"], result["error"], result["error"])
elif year == 2018:
    for etaKey, values in sorted(results[string1][string2].iteritems()):
        for ptKey, result in sorted(values.iteritems()):
            print "%s\t%s\t\t%.12f\t%.12f\t%.12f" % (etaKey, ptKey, result["value"], result["stat"], result["stat"])


print ("")
