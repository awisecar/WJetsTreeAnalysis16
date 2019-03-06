#! /usr/bin/env python2
import os
import sys
#from ROOT import *
import random
import time
import datetime

cwd = os.getcwd()
print 'Current working directory: ' + cwd + '\n'

os.system('root -b -q wjets_compileCode.cc')

####################################################################

print '\nCode finished compiling, writing Condor submit script!'

dateTo = datetime.datetime.now().strftime("%Y_%m_%d_%H%M%S")
mtmpdir = 'wjetsCondor_' + dateTo
os.system('mkdir -p ' + mtmpdir)
# os.system('cd '+mtmpdir)

#first make submit script, which needs a premade directory
submit = 'universe = vanilla\n'
# submit += 'executable = script.sh\n' #give executable argument on command line
submit += 'arguments = "$(argument)"\n'
submit += 'output = '+mtmpdir+'/wjetsSub_$(ClusterId)_$(ProcId).out\n'
submit += 'error = '+mtmpdir+'/wjetsSub_$(ClusterId)_$(ProcId).err\n'
submit += 'log = '+mtmpdir+'/wjetsSub_$(ClusterId)_$(ProcId).log\n\n'
submit += '+JobFlavour = "testmatch"\n\n'
submit += 'queue argument in 1'

submitName = mtmpdir+'_Submit.sub'
subScript = open(submitName,'w')
subScript.write(submit+'\n')
subScript.close()
os.system('chmod 755 '+submitName)


####################################################################

print 'Submit script finished, writing individual job scripts!'
cmsswdir = '/afs/cern.ch/user/a/awisecar/WJetsTreeAnalysis16/CMSSW_5_3_20/src'
# os.system('cd '+mtmpdir)

##doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 30, 41, 42, 51, 52, 53, 54, 61, 62, 63] #everything
#doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 30, 42, 51, 52, 53, 54] #important files
##doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19] #Data
## doWhat = [21, 22, 23, 24, 30] #Background
## doWhat = [41, 42, 51, 52, 53, 54, 61, 62, 63] #W+jets MC
#doWhat = [21, 22, 23, 24, 30, 41, 42, 51, 52, 53, 54, 61, 62, 63] #BG + W+jets MC for syst. uncert.'s
#
#doQCD = [0, 1, 2, 3] #signal + 3 control regions for QCD BG
##doQCD = [0]
#
#doSysRunning = [0] #nominal
##doSysRunning = [2] #JES uncertainties
##doSysRunning = [3, 4, 5, 6] #other uncertanties

#doWhat = [52]
#doQCD = [3]
#doSysRunning = [0]

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
doWhat = [51, 52, 53, 54] #W+jets MC
doQCD = [0]
doSysRunning = [4]

########### LepSF Syst
#doWhat = [21, 22, 23, 24, 30, 51, 52, 53, 54] #BG + W+jets MC for syst. uncert.'s
#doQCD = [0]
#doSysRunning = [5]
##############################

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
			
			com = 'root -b -q runDYJets.cc\(' + str(what) + ',' + str(QCD) + ',' + str(sys) + '\) 2>&1\n'

			jobName = mtmpdir+'/job_' + 'do' + str(what) + '_QCD' + str(QCD) + '_Sys' + str(sys) + '.sh'
			jobScript = open(jobName,'w')
			jobScript.write(job+'\n\n')
			jobScript.write(com)
			jobScript.close()
			os.system('chmod 755 '+jobName)

			command = 'condor_submit '+submitName+' executable='+jobName
			print '\nGoing to submit command: '+command
			print '-------> doWhat='+str(what)+', doQCD='+str(QCD)+', doSysRunning='+str(sys)

			os.system(command)
			os.system('sleep 1')


moveFile = 'mv '+submitName+' '+mtmpdir+'/'+submitName
print '\nCleaning up...'
print moveFile
os.system(moveFile)
print 'Finished!'

