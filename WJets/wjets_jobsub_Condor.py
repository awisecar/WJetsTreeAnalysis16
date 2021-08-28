#! /usr/bin/env python2
import os
import sys
import random
import time
import datetime

cwd = os.getcwd()
print '\nCurrent working directory: ' + cwd

####################################################################

print '\nWriting Condor submit script!'

dateTo = datetime.datetime.now().strftime("%Y_%m_%d_%H%M%S")
mtmpdir = 'wjetsCondor_' + dateTo
os.system('mkdir -p ' + mtmpdir)

#first make submit script, which needs a premade directory
#by default job gets one slot of a CPU core, 2Gb of memory and 20Gb of disk space
#can change this with "request_cpus" line below
submit  = 'universe = vanilla\n'
# submit += 'executable = script.sh\n' #give executable argument on command line
submit += 'arguments = "$(argument)"\n'
submit += 'output = '+mtmpdir+'/wjetsSub_$(ClusterId)_$(ProcId).out\n'
submit += 'error = '+mtmpdir+'/wjetsSub_$(ClusterId)_$(ProcId).err\n'
submit += 'log = '+mtmpdir+'/wjetsSub_$(ClusterId)_$(ProcId).log\n'
# submit += 'request_cpus = 2\n\n'

# submit += '+JobFlavour = "testmatch"\n\n' #testmatch is 3d queue
# submit += '+JobFlavour = "tomorrow"\n\n' #tomorrow is 1d queue
# submit += '+MaxRuntime = 43200\n\n' # set for 12h (12h = 43200s)
# submit += '+MaxRuntime = 36000\n\n' # set for 10h (10h = 36000s)
# submit += '+MaxRuntime = 32400\n\n' # set for 9h (9h = 32400s)
submit += '+JobFlavour = "workday"\n\n' #workday is 8h queue
# submit += '+JobFlavour = "espresso"\n\n' #espresso is 20min queue

submit += 'queue argument in 1'

submitName = mtmpdir+'_Submit.sub'
subScript = open(submitName,'w')
subScript.write(submit+'\n')
subScript.close()
os.system('chmod 755 '+submitName)

####################################################################

print 'Submit script finished, writing individual job scripts!'
### need to be able to grab cmsswdir automatically
cmsswdir = '/afs/cern.ch/user/a/awisecar/WJetsTreeAnalysis16_lxplus7/CMSSW_7_6_0/src'
# os.system('cd '+mtmpdir)

##############################
### ----- 2016 -----

#doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 221, 222, 223, 23, 24, 25, 26, 27, 30, 41, 42, 51, 52, 53, 54] # full set of files, w+jets pT-binned
# doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 221, 222, 223, 23, 24, 25, 26, 27, 30, 42, 61, 62, 63]   # full set of files, w+jets jet-binned
doWhat = [41]
#doWhat = [51, 52, 53, 54, 61, 62, 63]
# doWhat = [61, 62, 63]

#doQCD = [0, 1, 2, 3] # all regions
doQCD = [0] # signal region
#doQCD = [1, 2, 3] # control regions

doSysRunning = [0]

years = [2016]

##############################
### ----- 2017 -----

# doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 211, 212, 213, 221, 222, 223, 23, 24, 25, 26, 27, 30, 41, 42, 51, 52, 53, 54] # full set of files, w+jets pT-binned
# doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 211, 212, 213, 221, 222, 223, 23, 24, 25, 26, 27, 30, 42, 61, 62, 63]     # full set of files, w+jets jet-binned

# doQCD = [0, 1, 2, 3] # signal and control regions
# doQCD = [0] # signal region
# doQCD = [1, 2, 3] # control regions

# doSysRunning = [0]

# years = [2017]

##############################
### ----- 2018 -----

#doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 211, 212, 213, 221, 222, 223, 23, 24, 25, 26, 27, 30, 41, 42, 51, 52, 53, 54] # full set of files, w+jets pT-binned
# doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 211, 212, 213, 221, 222, 223, 23, 24, 25, 26, 27, 30, 42, 61, 62, 63]     # full set of files, w+jets jet-binned
#doWhat = [15]

#doQCD = [0, 1, 2, 3] #signal and control regions
#doQCD = [0] #signal region

#doSysRunning = [0]

#years = [2018]

##############################
### ----- Systematics -----

########## Pileup (PU)
#doWhat = [21, 221, 222, 223, 23, 24, 25, 26, 27, 30, 51, 52, 53, 54]            #Background & Signal, 2016
#doWhat = [211, 212, 213, 221, 222, 223, 23, 24, 25, 26, 27, 30, 51, 52, 53, 54] #Background & Signal, 2017+2018
#doQCD = [0]
#doSysRunning = [1]

########## Jet Energy Scale (JES)
#doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19] #Data
#doQCD = [0]
#doSysRunning = [2]

########## BG Cross Sections (XSEC)
#doWhat = [21, 221, 222, 223, 23, 24, 25, 26, 27, 30]            #Background, 2016
#doWhat = [211, 212, 213, 221, 222, 223, 23, 24, 25, 26, 27, 30] #Background, 2017+2018
#doQCD = [0]
#doSysRunning = [3]

########## Jet Energy Resolution (JER)
#doWhat = [51, 52, 53, 54] #W+jets MC
#doQCD = [0]
#doSysRunning = [4]

########### Lepton Eff. Scale Factors (LepSF)
#doWhat = [21, 221, 222, 223, 23, 24, 25, 26, 27, 30, 51, 52, 53, 54]            #Background & Signal, 2016
#doWhat = [211, 212, 213, 221, 222, 223, 23, 24, 25, 26, 27, 30, 51, 52, 53, 54] #Background & Signal, 2017+2018
#doQCD = [0]
#doSysRunning = [5]

########### B-Tagging Eff. Scale Factors (BTagSF) <-- NOT USED FOR NOW
#doWhat = [21, 221, 222, 223, 23, 24, 30, 51, 52, 53, 54]            #BG + W+jets MC for syst. uncert.'s, 2016
#doWhat = [211, 212, 213, 221, 222, 223, 23, 24, 30, 51, 52, 53, 54] #BG + W+jets MC for syst. uncert.'s, 2017+2018
#doQCD = [0]
#doSysRunning = [6]

########### L1 Prefiring Effect (L1Prefire)
#doWhat = [21, 221, 222, 223, 23, 24, 25, 26, 27, 30, 51, 52, 53, 54]            #Background & Signal, 2016
#doWhat = [211, 212, 213, 221, 222, 223, 23, 24, 25, 26, 27, 30, 51, 52, 53, 54] #Background & Signal, 2017+2018
#doQCD = [0]
#doSysRunning = [11]

##############################
### ----- ttbar SFs -----
# Remember to turn doBJets to 2!
# We do not run QCD BG for this control region study

# doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 211, 212, 213, 221, 222, 223, 23, 24, 30, 61, 62, 63]
# doQCD = [0]
# doSysRunning = [0]
# years = [2017]

##############################

for year in years:
	for what in doWhat:
		for QCD in doQCD:
			for sys in doSysRunning:
				
				job = '#!/bin/bash\n'
				job += 'export INPUT_ARG=$1\n\n'
				job += 'printf "Hello %s, your job is on host: $(hostname)'+str(r'\n')+'" "$USER"'+'\n'
				job += 'printf "Directory: $(pwd)'+str(r'\n\n')+'"\n'
				job += 'printf "Grabbing cmsenv and changing to WJetsTreeAnalysis workspace:'+str(r'\n')+'"\n'
				job += 'export CMSSW_PROJECT_SRC='+cmsswdir+'\n'
				job += 'cd $CMSSW_PROJECT_SRC\n'
				job += 'eval `scramv1 runtime -sh`\n'
				job += 'cd '+cwd+'\n'
				job += 'printf "$(pwd)'+str(r'\n')+'"\n'
				job += 'printf "CMSSW version: %s'+str(r'\n\n')+'" "$CMSSW_VERSION"\n\n'
				job += 'printf "Running code! ------------------------------- '+str(r'\n\n')+'"'
				
				com = 'root -l -b -q runEvtSelection.cc\(' + str(what) + ',' + str(QCD) + ',' + str(sys) + ',' + str(year) + '\) 2>&1\n'

				jobName = mtmpdir+'/job_' + str(year) + '_do' + str(what) + '_QCD' + str(QCD) + '_Sys' + str(sys) + '.sh'
				jobScript = open(jobName,'w')
				jobScript.write(job+'\n\n')
				jobScript.write(com)
				jobScript.close()
				os.system('chmod 755 '+jobName)

				command = 'condor_submit '+submitName+' executable='+jobName
				print '\nGoing to submit command: '+command
				print '-------> doWhat='+str(what)+', doQCD='+str(QCD)+', doSysRunning='+str(sys)+', year='+str(year)

				os.system(command)
				os.system('sleep 1')


moveFile = 'mv '+submitName+' '+mtmpdir+'/'+submitName
print '\nCleaning up...'
print moveFile
os.system(moveFile)

print '\nFinished!\n'
