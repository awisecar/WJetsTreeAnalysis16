import os
import sys
#from ROOT import *
import random
import time
import datetime

cwd = os.getcwd()
print 'Current working directory: ' + cwd + '\n'

os.system('root -b -q wjets_compileCode.cc')

dateTo = datetime.datetime.now().strftime("%Y_%m_%d_%H%M%S")
mtmpdir = 'wjetsjobs_' + dateTo
os.system('mkdir ' + mtmpdir)

doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 1, 3, 41, 42, 5, 6]
#doWhat = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
#doWhat = [1, 3, 41, 42, 5, 6]
#doWhat = [12, 18]
doQCD = [0, 1, 2, 3]
#doQCD = [0]
doSysRunning = [0]
#doSysRunning = [2]

cmsswdir = '/afs/cern.ch/user/a/awisecar/WJetsTreeAnalysis16/CMSSW_5_3_20/src'

print '\nCode finished compiling, beginning job submission:'

for what in doWhat:
	for QCD in doQCD:
		for sys in doSysRunning:

			tjobname_out = mtmpdir+'/job_' + 'do' + str(what) + '_QCD' + str(QCD) + '_Sys' + str(sys) + '.out'
			#tjobname_err = mtmpdir+'/job_' + 'do' + str(what) + '_QCD' + str(QCD) + '.err'
			tjobname = mtmpdir+'/job_' + 'do' + str(what) + '_QCD' + str(QCD) + '_Sys' + str(sys) + '.sh'

			job = '#!/bin/bash\n'
     	 	        #job += 'cd /afs/cern.ch/user/a/awisecar/WJetsTreeAnalysis16/CMSSW_5_3_20/src/WJetsTreeAnalysis16/WJets\n'
                        job += 'cd ' + cmsswdir + ' \n'
			job += 'eval `scramv1 runtime -sh` \n'
                        job += 'cd ' + cwd + ' \n'
		
			com = 'root -b -q runDYJets.cc\(' + str(what) + ',' + str(QCD) + ',' + str(sys) + '\) 2>&1'
			print '... going to submit run command ==> ', com
			com +='\n\n'

			ajob = str(job)
			tjob = open(tjobname,'w')
			tjob.write(ajob+'\n\n')
			tjob.write(com)
			tjob.close()
			os.system('chmod 755 '+tjobname)

        	        print '.out filename ==>', tjobname_out
          	        #print '.err filename ==>', tjobname_err

			#bsub = 'bsub -q 1nw -o ' +tjobname_out+ ' -e ' +tjobname_err+ ' -J ' +  tjobname + ' < ' + tjobname + ' '
      	                bsub = 'bsub -R "pool>30000" -q 2nd -o ' +tjobname_out+ ' -J ' +  tjobname + ' < ' + tjobname + ' '
			print bsub, '\n'
			os.system(bsub)
			os.system('sleep 1')





