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
submit = 'universe = vanilla\n'
# submit += 'executable = script.sh\n' #give executable argument on command line
submit += 'arguments = "$(argument)"\n'
submit += 'output = '+mtmpdir+'/wjetsSub_$(ClusterId)_$(ProcId).out\n'
submit += 'error = '+mtmpdir+'/wjetsSub_$(ClusterId)_$(ProcId).err\n'
submit += 'log = '+mtmpdir+'/wjetsSub_$(ClusterId)_$(ProcId).log\n\n'

##submit += '+JobFlavour = "testmatch"\n\n' #testmatch is 3d queue
# submit += '+JobFlavour = "tomorrow"\n\n' #tomorrow is 1d queue
# submit += '+JobFlavour = "workday"\n\n' #workday is 8h queue
submit += '+MaxRuntime = 43200\n\n' # set for 12h (12h = 43200s)
# submit += '+MaxRuntime = 36000\n\n' # set for 10h (10h = 36000s)
#submit += '+MaxRuntime = 32400\n\n' # set for 9h (9h = 32400s)

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

### Running full 2016 data/MC, incl. QCD BG
doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 30, 42, 51, 52, 53, 54, 61, 62, 63] # full set of files
#doWhat = [21]
#doWhat = [11]

doQCD = [0, 1, 2, 3]
# doQCD = [1,2,3]

doSysRunning = [0]

years = [2016]

##############################

# Running full 2017 data/MC, incl. QCD BG
#doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 211, 212, 213, 22, 23, 24, 25, 26, 27, 30, 61, 62, 63] # full set of files, w+jets jet-binned
#doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 211, 212, 213, 22, 23, 24, 25, 26, 27, 30, 52, 53, 54] # full set of files, w+jets pT-binned
# doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19] # only data
# doWhat = [211, 212, 213, 22, 23, 24, 25, 26, 27, 30, 61, 62, 63] # only MC, w+jets jet-binned
#doWhat = [13, 14, 23, 27, 30, 63]
#doWhat = [23, 27, 63]
# doWhat = [211, 212, 213, 22, 23, 24, 30, 61, 62, 63] # just MC
#doWhat = [61, 62, 63] # just W+jets MC
# doWhat = [213]

#doQCD = [0, 1, 2, 3] # signal and control regions
# doQCD = [0] # signal region
# doQCD = [1, 2, 3] # control regions
# doQCD = [2,3]

# doSysRunning = [0]

# years = [2017]

##############################

## Running full 2018 data/MC, incl. QCD BG
# doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 211, 212, 213, 22, 23, 24, 25, 26, 27, 30, 61, 62, 63] # full set of files
# # doWhat = [211, 212, 213]
# # doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 30, 61, 62, 63] 
# # doWhat = [213]
# doQCD = [0, 1, 2, 3] #signal and control regions
# # doQCD = [0, 1]
# # doQCD = [0] #signal region
# doSysRunning = [0]
# years = [2018]

##############################
## Systematics ---

########## PU Syst
#doWhat = [21, 22, 23, 24, 30, 51, 52, 53, 54] #BG + W+jets MC for syst. uncert.'s
#doQCD = [0]
#doSysRunning = [1]

########## JES Syst
#doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19] #Data
#doQCD = [0]
#doSysRunning = [2]

########## XSec Syst
#doWhat = [21, 22, 23, 24, 30] #Background
#doQCD = [0]
#doSysRunning = [3]

########## JER Syst
#doWhat = [51, 52, 53, 54] #W+jets MC
#doQCD = [0]
#doSysRunning = [4]

########### LepSF Syst
#doWhat = [21, 22, 23, 24, 30, 51, 52, 53, 54] #BG + W+jets MC for syst. uncert.'s
#doQCD = [0]
#doSysRunning = [5]

########### BTagSF Syst
#doWhat = [21, 22, 23, 24, 30, 51, 52, 53, 54] #BG + W+jets MC for syst. uncert.'s
#doQCD = [0]
#doSysRunning = [6]

##############################
## Migrations Study
#doWhat = [41, 51, 52, 53, 54]
#doQCD = [0]
#doSysRunning = [0]
##############################

##############################
## ttbar SFs (remember to turn doBJets to 2)
## we do not run QCD BG for this control region study
# doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 211, 212, 213, 22, 23, 24, 30, 61, 62, 63]
# doWhat = [61, 62, 63]
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

