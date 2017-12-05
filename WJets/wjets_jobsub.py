
import os
import sys
#from ROOT import *
import random
import time
import datetime

os.system('root -b -q wjets_compileCode.cc')

dateTo = datetime.datetime.now().strftime("%Y_%m_%d_%H%M%S")
mtmpdir = 'wjetsjobs_' + dateTo
os.system('mkdir ' + mtmpdir)

#doWhat = [0,1,3,4]
#The 101...110 is for the data
doWhat = [101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 1, 3, 4]
#doWhat = [110]
#doQCD = [0]
#doQCD = [1,2,3]
doQCD = [0,1,2,3]


for what in doWhat:
	for QCD in doQCD:

		tjobname_out = mtmpdir+'/job_' + 'do' + str(what) + '_QCD' + str(QCD) + '.out'
		tjobname_err = mtmpdir+'/job_' + 'do' + str(what) + '_QCD' + str(QCD) + '.err'
		tjobname = mtmpdir+'/job_' + 'do' + str(what) + '_QCD' + str(QCD) + '.sh'

		job = '#!/bin/bash\n'
		job += 'cd $CMSSW_BASE/src/WJetsTreeAnalysis16/WJets\n'
		job += 'eval `scramv1 runtime -sh`'
		#job += 'cd $CMSSW_BASE/src\n'
		#job += 'cmsenv\n'
		#job += 'mycwd=`pwd`\n'
		#job += 'cd $mycwd\n'
		
		com = 'root -b -q runDYJets.cc\(' + str(what) + ',' + str(QCD) + '\)'
		print '... going to submit run command ==> ', com
		com +='\n\n'

		ajob = str(job)
		tjob = open(tjobname,'w')
		tjob.write(ajob+'\n\n')
		tjob.write(com)
		tjob.close()
		os.system('chmod 755 '+tjobname)

		bsub = 'bsub -q 2nw -o ' +tjobname_out+ ' -e ' +tjobname_err+ ' -J ' +  tjobname + ' < ' + tjobname + ' '
		print bsub, '\n'
		os.system(bsub)
		os.system('sleep 1')





