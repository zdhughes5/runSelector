#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 13:37:54 2021

@author: zdhughes
"""

#         MASTER FILE                            CHILD FILE
#############################          #############################
#[General]                  #          # [StageConfig]             #
#>Options                   #          # >settings                 #
#                           #          #                           #
#[Sections]                 #          #                           #
#>selections                #          #                           #
#                           #          #                           #
#[Selection]                #          #                           #
#>stages                    #          #                           #
#                           #          #                           #
#[SelectionConfig]          #          #                           #
#>parameters                #          #                           #
#                           #          #                           #
#                           #          #                           #
#############################          #############################

import general_functions_ed as GFF
from configparser import ConfigParser, ExtendedInterpolation
import subprocess
import sys
import code
import os
import numpy as np
from collections import OrderedDict as od
import xml.etree.ElementTree as ET
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table


#This class object should hold all the info for a "top" stage. I.e. "Main"+"MainConfig"+the main child config.
#It gets passed the masterConfig OBJECT and a name of a section within that config object.
#It pulls out that section and the section+"config" section and packages it as a single object.
class runlistMember():
	
	"""This holds all the information related to a particular run as defined in the master configuration file."""
	
	def __init__(self, masterConfig = None, selectionName = None, c = None):
		
		"""Loads in the two related sections from the master configuration file and loads in the child config file
		that holds all the parameters."""
		
		self.masterConfig = masterConfig if masterConfig is not None else sys.exit('No master config file provided.')
		self.selectionName = selectionName if selectionName is not None else sys.exit('No top section to load provided.')
		self.c = c if c is not None else GFF.colors()
		self.selectionConfigName = self.selectionName+'Config'

		#From master:
		#[General] options -> general
		#The selection options -> stages
		#the selection config options -> selection config
		self.general = od(self.masterConfig.items('General'))
		self.stages = od(self.masterConfig.items(self.selectionName))
		self.parameters = od(self.masterConfig.items(self.selectionConfigName))

		#code.interact(local=locals())
		#sys.exit('Code Break!')	
		
		#Cast the above into the correct types.
		for key in self.general:
			self.general[key] = GFF.tryeval(self.general[key])
		for stage in self.stages:
			self.stages[stage] = GFF.tryeval(self.stages[stage])	
		for parameter in self.parameters:
			self.parameters[parameter] = GFF.tryeval(self.parameters[parameter])

		#Read in the child config, as defined in the selection config. DO NOT cast them as types, since they will be written out to files.
		with open(os.path.expandvars(self.parameters['anasumConfig']), 'r') as f:
			self.childConfig = [x.strip('\n') for x in f]
			
	def printInfo(self):
		
		print(self.c.gray('++++++++++++++++++++++++++++++'))
		
		print('Printing section '+self.c.green(self.selectionName))
		for stage in self.stages:
			if self.stages[stage] == True:
				print(self.c.lightgreen(stage)+' : '+self.c.lightblue(self.stages[stage]))
			else:
				print(self.c.lightgreen(stage)+' : '+self.c.lightred(self.stages[stage]))

		print(self.c.orange('------------------------------'))
		
		print('Printing config section '+self.c.green(self.selectionConfigName))
		for parameter in self.parameters:
			print(self.c.lightgreen(parameter)+' : '+self.c.lightblue(self.parameters[parameter]))
		
		print(self.c.orange('------------------------------'))

		print('Printing active anasum parameters '+self.c.yellow(self.parameters['anasumConfig']))
		for line in self.childConfig:
			if line and line[0] == '*':
				leading = ' '.join(line.split()[0:2]) + ' '
				trailing = ' '.join(line.split()[2:])		
				print(self.c.lightgreen(leading) + ' : ' + self.c.lightblue(trailing))
				if leading + trailing != line:
					print(self.c.red('Warning did not reproduce anasum config line correctly!'))
					print('Leading:')
					print(leading)
					print('Trailing:')
					print(trailing)
					print('Line:')
					print(line)
			else:
				continue
				print(self.c.lightgreen(line))
				
		print(self.c.gray('++++++++++++++++++++++++++++++'))
				
	def printAnaConfig(self):
		
		print(self.c.gray('++++++++++++++++++++++++++++++'))
		
		print('Printing anasum config file ' + self.c.yellow(self.parameters['anasumConfig']))
		for line in self.childConfig:
			if line and line[0] == '*':
				leading = ' '.join(line.split()[0:2]) + ' '
				trailing = ' '.join(line.split()[2:])		
				print(self.c.lightgreen(leading) + ' : ' + self.c.lightblue(trailing))
				if leading + trailing != line:
					print(self.c.red('Warning did not reproduce anasum config line correctly!'))
					print('Leading:')
					print(leading)
					print('Trailing:')
					print(trailing)
					print('Line:')
					print(line)
			else:
				print(self.c.lightgreen(line))
				
		print(self.c.gray('++++++++++++++++++++++++++++++'))
		
	def getChild(self):
		return self.childConfig
	

class EDRunlistMember(runlistMember):
	
	def __init__(self, masterConfig = None, selectionName = None, c = None):
		runlistMember.__init__(self, masterConfig, selectionName, c)
		with open(self.general['runlistFile'], 'r') as f:
			self.runs = [x.strip('\n') for x in f]
		self.general['saveDir'] = os.path.expandvars(self.general['saveDir'])
		self.parameters['anasumConfig'] = os.path.expandvars(self.parameters['anasumConfig'])
			
	def __iter__(self):
		return iter(self.times)
			

class ecumenical():
	
	"""This object's primary purpose is to be a container for runlistMember objects."""
	
	def __init__(self, root_python_program = None, masterConfigFilename = None, c = None):
		
		"""Given a path to the master config file, Pretty much everything useful is done at instantiation."""
		
		self.root_python_program = root_python_program if root_python_program is not None else '##DEFAULT##'
		self.masterConfigFilename = masterConfigFilename if masterConfigFilename is not None else sys.exit('Error. No master Configuration file!')
		self.c = c if c is not None else GFF.colors()
	
		self.masterConfig = ConfigParser(interpolation=ExtendedInterpolation(),inline_comment_prefixes=('#'))
		self.masterConfig.optionxform = str
		self.masterConfig.read(masterConfigFilename)
		
		print('Loading in the master configuration file...')
		
		self.general = od(self.masterConfig['General'])
		self.sections = od(self.masterConfig['Sections'])

		for key in self.general:
			self.general[key] = GFF.tryeval(self.general[key])		
		for selection in self.sections:
			self.sections[selection] = GFF.tryeval(self.sections[selection])		
		
		if self.general['ansiColors'] == True:
			self.c.enableColors()
			self.c.confirmColorsDonger()

		self.members = od({})
		for selection in self.sections:
			if self.sections[selection] == True:
				self.members[selection] = EDRunlistMember(self.masterConfig, selection, self.c)
				
	def __getitem__(self, key):
		return self.members[key]

	def to_list(self):
		return list(self.members.values())

class collectModels():
	pass	
	