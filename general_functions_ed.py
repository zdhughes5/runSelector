#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 13:26:45 2021

@author: zdhughes
"""

from time import gmtime, strftime, sleep, time
import sys
import ast
import subprocess
import code
import os

# v483b
# =============================================================================
# eventDisplayJobPath = '/data/programs/collaborations/veritas/EventDisplay/v483b/scripts/VTS/ANALYSIS.evndisp.sh'
# lookupJobPath = '/data/programs/collaborations/veritas/EventDisplay/v483b/scripts/VTS/ANALYSIS.mscw_energy.sh'
# anasumJobPath = '/data/programs/collaborations/veritas/EventDisplay/v483b/scripts/VTS/ANALYSIS.anasum_parallel_from_runlist.sh'
# combineJobPath = '/data/programs/collaborations/veritas/EventDisplay/v483b/scripts/VTS/ANALYSIS.anasum_combine.sh'
# =============================================================================

eventDisplayJobPath = '/data/programs/collaborations/veritas/EventDisplay/v485a/Eventdisplay_AnalysisScripts_VTS/scripts/ANALYSIS.evndisp.sh'
lookupJobPath = '/data/programs/collaborations/veritas/EventDisplay/v485a/Eventdisplay_AnalysisScripts_VTS/scripts/ANALYSIS.mscw_energy.sh'
anasumJobPath = '/data/programs/collaborations/veritas/EventDisplay/v485a/Eventdisplay_AnalysisScripts_VTS/scripts/ANALYSIS.anasum_parallel_from_runlist.sh'
combineJobPath = '/data/programs/collaborations/veritas/EventDisplay/v485a/Eventdisplay_AnalysisScripts_VTS/scripts/ANALYSIS.anasum_combine.sh'


def getCondorEDParameters(directory, sectionName, run):
	
	"""Given a gttool section and location this function will return parameters central to creating a fermi condor script."""
	
	if sectionName == 'doEventDisplay':
		job = directory+'condor_eventDisplayJob'
		#executable = '$EVNDISPSYS/scripts/VTS/ANALYSIS.evndisp.sh'
		executable = eventDisplayJobPath
		textDir = directory+'eventDisplayText'
		entries = {
			'output' : textDir+'/eventDisplay_'+run+'.out',
			'log' : textDir+'/eventDisplay_'+run+'.log',
			'error' : textDir+'/eventDisplay_'+run+'.err'
			}

	if sectionName == 'doLookup':
		job = directory+'condor_lookupJob'
		#executable = '$EVNDISPSYS/scripts/VTS/ANALYSIS.mscw_energy.sh'
		executable = lookupJobPath
		textDir = directory+'lookupText'
		entries = {
			'output' : textDir+'/lookup_'+run+'.out',
			'log' : textDir+'/lookup_'+run+'.log',
			'error' : textDir+'/lookup_'+run+'.err'
			}

	if sectionName == 'doAnasum':
		job = directory+'condor_anasumJob'
		#executable = '$EVNDISPSYS/scripts/VTS/ANALYSIS.anasum_parallel_from_runlist.sh'
		executable = anasumJobPath
		textDir = directory+'anasumText'
		entries = {
			'output' : textDir+'/anasum.out',
			'log' : textDir+'/anasum.log',
			'error' : textDir+'/anasum.err'
			}
		
	if sectionName == 'doCombine':
		job = directory+'condor_combineJob'
		executable = combineJobPath
		textDir = directory+'combineText'
		entries = {
			'output' : textDir+'/combine.out',
			'log' : textDir+'/combine.log',
			'error' : textDir+'/combine.err'
			}
		

	return job, executable, textDir, entries

class tee():
	
	"tee implemented in python"
	
	def __init__(self, fd1, fd2):
		self.fd1 = fd1
		self.fd2 = fd2
		
	def __del__(self):
		
		stdStreams = [0, 1, 2]
		fd1_no = self.fd1.fileno()
		fd2_no = self.fd2.fileno()
		
		if fd1_no not in stdStreams:
			self.fd1.close()
		if fd2_no not in stdStreams:
			self.fd2.close()
	
	def write(self, text):
		self.fd1.write(text)
		self.fd2.write(text)

	def flush(self):
		self.fd1.flush()
		self.fd2.flush()


class colors():
	
	"""This class will color console text output using ANSI codes."""

	RED = ''
	ORANGE = ''
	YELLOW = ''
	GREEN = ''
	BLUE = ''
	INDIGO = ''
	VIOLET = ''
	PINK = ''                                                                        
	BLACK = ''
	CYAN = ''
	PURPLE = ''
	BROWN = ''
	GRAY = ''
	DARKGRAY = ''
	LIGHTBLUE = ''
	LIGHTGREEN = ''
	LIGHTCYAN = ''
	LIGHTRED = ''
	LIGHTPURPLE = ''
	WHITE = ''
	BOLD = ''
	UNDERLINE = ''
	ENDC = ''
    
	def enableColors(self):

		self.RED = '\033[0;31m'
		self.ORANGE = '\033[38;5;166m'
		self.YELLOW = '\033[1;33m'
		self.GREEN = '\033[0;32m'
		self.BLUE = '\033[0;34m'
		self.INDIGO = '\033[38;5;53m'
		self.VIOLET = '\033[38;5;163m'
		self.PINK =  '\033[38;5;205m'
		self.BLACK = '\033[0;30m'
		self.CYAN = '\033[0;36m'
		self.PURPLE = '\033[0;35m'
		self.BROWN = '\033[0;33m'
		self.GRAY = '\033[0;37m'
		self.DARKGRAY = '\033[1;30m'
		self.LIGHTBLUE = '\033[1;34m'
		self.LIGHTGREEN = '\033[1;32m'
		self.LIGHTCYAN = '\033[1;36m'
		self.LIGHTRED = '\033[1;31m'
		self.LIGHTPURPLE = '\033[1;35m'
		self.WHITE = '\033[1;37m'
		self.BOLD = '\033[1m'
		self.UNDERLINE = '\033[4m'
		self.ENDC = '\033[0m'

        
	def disableColors(self):
        
		self.RED = ''        
		self.ORANGE = ''
		self.YELLOW = ''
		self.GREEN = ''
		self.BLUE = ''
		self.INDIGO = ''
		self.VIOLET = ''
		self.PINK = ''
		self.BLACK = ''
		self.CYAN = ''
		self.PURPLE = ''
		self.BROWN = ''
		self.GRAY = ''
		self.DARKGRAY = ''
		self.LIGHTBLUE = ''
		self.LIGHTGREEN = ''
		self.LIGHTCYAN = ''
		self.LIGHTRED = ''
		self.LIGHTPURPLE = ''
		self.WHITE = ''
		self.BOLD = ''
		self.UNDERLINE = ''
		self.ENDC = ''

	def getState(self):
		if self.ENDC:
			return True
		elif not self.ENDC:
			return False
		else:
			return -1

	def flipState(self):
		if self.getState():
			self.disableColors()
		elif not self.getState():
			self.enableColors()
		else:
			sys.exit("Can't flip ANSI state, exiting.")

	def confirmColors(self):
		if self.getState() == True:
			print('Colors are '+self.red('e')+self.orange('n')+self.yellow('a')+self.green('b')+self.blue('l')+self.indigo('e')+self.violet('d'))
		elif self.getState() == False:
			print('Colors are off!')
		elif self.getState() == -1:
			print('Error: Can\'t get color state.')

	def confirmColorsDonger(self):
		if self.getState() == True:
			print('Colors are '+self.pink('(ﾉ')+self.lightblue('◕')+self.pink('ヮ')+self.lightblue('◕')+self.pink('ﾉ')+self.red('☆')+self.orange('.')+self.yellow('*')+self.green(':')+self.blue('･ﾟ')+self.indigo('✧')+self.violet(' enabled!'))    
		elif self.getState() == False:
			print('Colors are off!')
		elif self.getState() == -1:
			print('Error: Can\'t get color state.')

	def orange(self, inString):
		inString = str(self.ORANGE+str(inString)+self.ENDC)
		return inString
	def indigo(self, inString):
		inString = str(self.INDIGO+str(inString)+self.ENDC)
		return inString
	def violet(self, inString):
		inString = str(self.VIOLET+str(inString)+self.ENDC)
		return inString
	def pink(self, inString):
		inString = str(self.PINK+str(inString)+self.ENDC)
		return inString
	def black(self, inString):
		inString = str(self.BLACK+str(inString)+self.ENDC)
		return inString
	def blue(self, inString):
		inString = str(self.BLUE+str(inString)+self.ENDC)
		return inString
	def green(self, inString):
		inString = str(self.GREEN+str(inString)+self.ENDC)
		return inString
	def cyan(self, inString):
		inString = str(self.CYAN+str(inString)+self.ENDC)
		return inString
	def red(self, inString):
		inString = str(self.RED+str(inString)+self.ENDC)
		return inString
	def purple(self, inString):
		inString = str(self.PURPLE+str(inString)+self.ENDC)
		return inString
	def brown(self, inString):
		inString = str(self.BROWN+str(inString)+self.ENDC)
		return inString
	def gray(self, inString):
		inString = str(self.GRAY+str(inString)+self.ENDC)
		return inString
	def darkgray(self, inString):
		inString = str(self.DARKGRAY+str(inString)+self.ENDC)
		return inString
	def lightblue(self, inString):
		inString = str(self.LIGHTBLUE+str(inString)+self.ENDC)
		return inString
	def lightgreen(self, inString):
		inString = str(self.LIGHTGREEN+str(inString)+self.ENDC)
		return inString
	def lightcyan(self, inString):
		inString = str(self.LIGHTCYAN+str(inString)+self.ENDC)
		return inString
	def lightred(self, inString):
		inString = str(self.LIGHTRED+str(inString)+self.ENDC)
		return inString
	def yellow(self, inString):
		inString = str(self.YELLOW+str(inString)+self.ENDC)
		return inString
	def white(self, inString):
		inString = str(self.WHITE+str(inString)+self.ENDC)
		return inString
	def bold(self, inString):
		inString = str(self.BOLD+str(inString)+self.ENDC)
		return inString
	def underline(self, inString):
		inString = str(self.UNDERLINE+str(inString)+self.ENDC)
		return inString

class condorHandler():
	
	"""Object used to open, write, and run a condor vanilla-universe submit file."""
	
	def __init__(self, filename=None, executable=None, title=None, subtitle=None, tertiary=None, c=None):
		
		"""Opens a write job at the filename location. Writes the initial header."""
		
		self.c = c if c is not None else colors()
		colorState = self.c.getState()
		self.filename = filename if filename is not None else sys.exit(c.red('Error: Job filename must be given.'))
		if executable == None: sys.exit(c.red('Error: Job executable must be given.'))
		self.condorJob = open(filename,'w')
		self.logs = []
	
		if colorState: c.disableColors()
		self.condorJob.write('######################################\n')
		self.condorJob.write('#\n')
		self.condorJob.write('# '+title+'\n')
		self.condorJob.write('# '+subtitle+'\n')
		self.condorJob.write('# '+tertiary+'\n')
		self.condorJob.write('# Generated at UTC '+strftime("%Y-%m-%d %H:%M:%S", gmtime())+'\n')
		self.condorJob.write('#\n')
		self.condorJob.write('######################################\n')
		self.condorJob.write('\n')
		self.condorJob.write('Universe   = vanilla\n')
		self.condorJob.write('Executable = '+executable+'\n')
		self.condorJob.write('\n')
		self.condorJob.flush()
		if colorState: c.enableColors()

	def addEntry(self, entries):
		
		"""Takes a dictionary and adds them as condor queue entries."""
		
		colorState = self.c.getState()
	
		if entries == None:
			print(self.c.red('Warning: no entries given in condor header creator. This is probably wrong!'))
	
		if colorState: self.c.disableColors()	
		if 'requirements' in entries:
			self.condorJob.write('requirements'+' = '+entries['requirements']+'\n')
			del entries['requirements']
		if 'output' in entries:
			self.condorJob.write('output'+' = '+entries['output']+'\n')
			del entries['output']
		if 'error' in entries:
			self.condorJob.write('error'+' = '+entries['error']+'\n')
			del entries['error']
		if 'log' in entries:
			self.condorJob.write('log'+' = '+entries['log']+'\n')
			self.logs.append(entries['log'])
			del entries['log']
		if 'arguments' in entries:
			self.condorJob.write('arguments'+' = '+entries['arguments']+'\n')
			del entries['arguments']
	
		if len(entries) > 0:
			for key in entries:
				self.condorJob.write(key+' = '+entries[key]+'\n')
				
		self.condorJob.write('getenv = True\n')
		self.condorJob.write('queue\n')
		self.condorJob.write('\n')
		self.condorJob.flush()
	
		if colorState: self.c.enableColors()
				
	def close(self):
		self.condorJob.close()
		
	def run(self, timeout = None, prepend = '', append = '', verbosity=3):
		submissionCommand = prepend + 'condor_submit ' + self.filename + append
		if verbosity >= 1: print('Running this command: ', submissionCommand)	
		if verbosity >= 2: print('Running with this timeout: ', timeout)
		if verbosity >= 3: print('Currently in directory: ', os.getcwd())
		if verbosity >= 3: print('Looking at these logs: ', self.logs)
		if timeout:
			self.starttime = time()
		thing = subprocess.call(submissionCommand, shell=True, executable='/bin/bash')
		#print(thing)
		while len(self.logs) > 0:
			for i, log in enumerate(self.logs):
				if subprocess.run('grep "return" '+log, shell=True).returncode == 0:
					del self.logs[i]
			sleep(1)
			if timeout:
				if (time() - self.starttime) > timeout:
					print(self.c.yellow('Timeout reached. Breaking.'))
					break
		else:
			if verbosity >= 3: print('Number of log files has dropped below 1: ', self.logs)
	
def convertSectionTypes(section):

	"""Converts a dictionary to their most natural types."""
	
	for i, key in enumerate(section):
		try:
			section[key] = ast.literal_eval(section[key])
		except (SyntaxError, ValueError) as e:
			pass
	return section

### Old, do not use. ###

def makeCondorHeader(job=None, executable=None, program=None, stage=None, condorFile=None, c=None):

	if c == None:
		c = colors()

	colorState = c.getState()

	if executable == None:
		sys.exit(c.red('Warning: no executable given in condor header creator.'))
	if job == None:
		sys.exit(c.red('Warning: no job file given in condor header creator.'))

	if colorState:
		c.disableColors()

	job.write('######################################\n')
	job.write('#\n')
	job.write('# '+program+'\n')
	job.write('# '+stage+'\n')
	job.write('# '+condorFile+'\n')
	job.write('# Generated at UTC '+strftime("%Y-%m-%d %H:%M:%S", gmtime())+'\n')
	job.write('#\n')
	job.write('######################################\n')
	job.write('\n')
	job.write('Universe   = vanilla\n')
	job.write('Executable = '+executable+'\n')
	job.write('\n')

	if colorState:
		c.enableColors()

def makeCondorBody(job=None, entries=None, c=None):

	if c == None:
		c = colors()

	colorState = c.getState()

	if job == None:
		sys.exit(c.red('Warning: no job file given in condor header creator.'))
	if entries == None:
		print(c.red('Warning: no entries given in condor header creator. This is probably wrong!'))

	if colorState:
		c.disableColors()

	if 'requirements' in entries:
		job.write('requirements'+' = '+entries['requirements'])
		del entries['requirements']
	if 'output' in entries:
		job.write('output'+' = '+entries['output'])
		del entries['output']
	if 'error' in entries:
		job.write('error'+' = '+entries['error'])
		del entries['error']
	if 'log' in entries:
		job.write('log'+' = '+entries['log'])
		del entries['log']
	if 'arguments' in entries:
		job.write('arguments'+' = '+entries['arguments'])
		del entries['arguments']

	if len(entries) > 0:
		for key in entries:
			job.write(key+' = '+entries[key])
	job.write('getenv = True\n')
	job.write('queue\n')
	job.write('\n')

	if colorState:
		c.enableColors()

def setCommandlineArgs(member, stage, run):

	lineToReturn = ''
	if stage == 'doEventDisplay':
		lineToReturn = member.parameters['directory']+'TMPDIR/'+run+'.lst'
	elif stage == 'doLookup':
		lineToReturn = member.parameters['directory']+'TMPDIR/'+run+'.lst '
		lineToReturn += member.general['saveDir']
	elif stage == 'doAnasum':
		lineToReturn = member.general['runlistFile'] + ' '
		lineToReturn += member.general['saveDir'] + member.parameters['saveNameStem'] + ' '
		lineToReturn += member.parameters['cutSet'] + ' '
		lineToReturn += member.parameters['backgroundModel'] + ' '
		lineToReturn += member.parameters['anasumConfig']
	elif stage == 'doCombine':
		lineToReturn = member.general['saveDir'] + member.parameters['saveNameStem'] + '/' + member.parameters['cutSet'] + '.anasum.dat '
		lineToReturn += member.general['saveDir'] + member.parameters['saveNameStem'] + ' '
		lineToReturn += member.parameters['saveNameStem'] + '.' + member.parameters['backgroundModel'] + '.root '
		lineToReturn += member.parameters['anasumConfig']
	
	return lineToReturn


def tryeval(val):
	try:
		val = ast.literal_eval(val)
	except:
		pass
	return val







