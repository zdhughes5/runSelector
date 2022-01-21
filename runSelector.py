#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 18:32:53 2022

@author: zdhughes
"""

import numpy as np
import ROOT
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cf 
from configparser import ConfigParser, ExtendedInterpolation
import subprocess
import matplotlib.cm as cm 
import math
from ast import literal_eval as le
import scipy.ndimage.interpolation as ni
from matplotlib.patches import Rectangle
import pandas as pd
from multiprocessing import Pool, Lock, Manager
import os.path
from astropy.time import Time
from collections import Counter
from matplotlib import rcParams
from scipy import signal
import scipy.interpolate as si



def distance(x, popt, xPoint, yPoint):

	return np.sqrt((x - xPoint)**2 + (cosineFit(x, *popt) - yPoint)**2)

def d2r(x):
	
	return (np.pi/180) * x

def cosineFit(x, A):
	
	return A * np.cos(d2r(x)) ** 0.9

def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

def binRate(timeArray, pool, interval=10):
	
	startTime = timeArray[0]
	endTime = timeArray[-1]
	duration = endTime - startTime
	
	nArrays = int(duration/interval) + 1
	
	subarrs = np.array_split(timeArray, nArrays, axis=0)

	results = [pool.apply_async(calculateRate, args=(subarr,)) for subarr in subarrs]
	output = [p.get() for p in results]
	binnedRate, binnedRateErr = list(zip(*output))
	
	return (np.array(binnedRate), np.array(binnedRateErr))

def getInterpLiveTime(binnedLiveTime, finalArrSize):
	
	zoomFactor = float(finalArrSize/binnedLiveTime.shape[0])
	zoomedLivetime = ni.zoom(binnedLiveTime, zoomFactor)
	return zoomedLivetime

def splineCorrection2(wf, finalArrSize=4480, s=300, nDivisions = 50):
	
	splitwfs = np.array_split(wf, nDivisions)
	inds = []
	offset_precusor = [0] + [x.shape[0] for x in splitwfs]
	offsets = np.cumsum(offset_precusor)
	
	for i, splitwf in enumerate(splitwfs):
		thisMed = np.median(splitwf)
		thisInd = np.argmin(np.abs(splitwf - thisMed)) + offsets[i]
		inds.append(thisInd)
		
	x = np.array(inds)
	y = np.array(wf[x])
	tck = si.splrep(x, y, s=s)
	yFit = si.splev(x, tck, der=0)
	zoomFactor = float(finalArrSize/len(x))
	yZoomed = ni.zoom(yFit, zoomFactor)

	return yZoomed

def calculateRate(times):
	
	#print('Calculating rate...')
	rate = np.nan
	rateErr = np.nan
	
	timesZ = times - times[0]
	
	time_diffs = np.diff(timesZ)
	
	# Too lazy to use np.digitize
	f, ax = plt.subplots()
	n, bins, patches = ax.hist(time_diffs, bins=1000, density=True)
	plt.close()
	
	# Center points
	halfWidth = (bins[-2] - bins[-1])/2
	bins = bins[0:-1] - halfWidth
	p0 = [150]
	xData = bins
	yData = n

	# Fit
	try:
		popt, pcov = cf(timeIntervalFunction, xData, yData, p0=p0)
		rate = popt[0]
		rateErr = np.sqrt(np.diagonal(pcov))[0]
	except Exception as e:
		print(e)

	return rate, rateErr	

def divide_chunks(l, n):
      
	# looping till length l
	for i in range(0, len(l), n): 
		yield l[i:i + n]

def generateSeasonDF():
	# Dates pulled from ED's VERITAS.Epochs.runparameter file.
	dates = list(divide_chunks(Time([
		'2006-11-09',
		'2007-11-24',
		'2008-11-12',
		'2009-11-02',
		'2010-11-20',
		'2011-11-10',
		'2012-10-29',
		'2013-11-17',
		'2014-11-08',
		'2015-11-26',
		'2016-11-15',
		'2017-12-06',
		'2018-10-25',
		'2019-11-13',
		'2020-11-05',
		'2007-06-03',
		'2008-05-18',
		'2009-06-05',
		'2010-05-26',
		'2011-05-15',
		'2012-05-04',
		'2013-05-22',
		'2014-05-12',
		'2015-05-31',
		'2016-05-20',
		'2017-05-10',
		'2018-04-30',
		'2019-05-18',
		'2020-05-07',
		'2021-04-27',
		'2007-06-04',
		'2008-05-19',
		'2009-06-06',
		'2010-05-27',
		'2011-05-16',
		'2012-05-05',
		'2013-05-23',
		'2014-05-13',
		'2015-06-01',
		'2016-05-21',
		'2017-05-11',
		'2018-05-01',
		'2019-05-19',
		'2020-05-08',
		'2021-04-28',
		'2007-11-23',
		'2008-11-11',
		'2009-11-01',
		'2010-11-19',
		'2011-11-09',
		'2012-10-28',
		'2013-11-16',
		'2014-11-07',
		'2015-11-25',
		'2016-11-14',
		'2017-12-05',
		'2018-10-24',
		'2019-11-12',
		'2020-11-04',
		'2021-11-04'
		]).sort(), 2))
	
	
	# Define the Winter and Summer seasons
	dateIndex = [str(2000 + int(x/2)) + 'S' if x % 2 == 0 else str(2000 + int(x/2)) + 'W' for x in range(13, 43)]
	seasonsDF = pd.DataFrame({'Dates' : dates}, index=dateIndex)
	return seasonsDF

def timeIntervalFunction(t, r):
	
	return r * np.exp(-r*t)

class mscwFile:
	
	''' A class to open and hold a root file. Can calculate the relevent event rate. '''
	
	def __init__(self, filename, verbosity = 0):
		
		self.filename = filename
		self.verbosity = verbosity
		self.fileOpen = False
		self.rate = None
		self.rateErr = None
		self.times = None
		self.elevation = None
		self.drop = None
		self.binnedLiveTime = None
		
	def load(self):
		
		self.rootFile = ROOT.TFile.Open(self.filename, "READ")
		try:
			if self.rootFile.GetPath(): self.fileOpen = True		
			if self.verbosity >= 2:
				print('Opened file ', self.filename)
		except ReferenceError:
			self.fileOpen = False
			if self.verbosity >= 1:
				print('Failed to open file at ', self.filename)

	
	def extractEventRate(self):
		if self.fileOpen:
			self.rate, self.rateErr, self.times, self.elevation, self.drop, self.binnedLiveTime = self.calculateEventRate(self.rootFile)
		else:
			if self.verbosity >= 1:
				print('File not open.')
				
	def getRateInfo(self):
		return self.rate, self.rateErr, self.times, self.elevation, self.drop, self.binnedLiveTime
				
	
	@staticmethod
	def calculateEventRate(rootFile):
		
		# Pre-allocate
		times = np.zeros(rootFile.data.GetEntriesFast())
		rate = None
		elevations = np.zeros(rootFile.data.GetEntriesFast())
		drop = False
		
		# Loop over events, get times.
		for event in range(rootFile.data.GetEntriesFast()):
			rootFile.data.GetEntry(event)
			times[event] = rootFile.data.Time
			elevations[event] = rootFile.data.TelElevation[0]
			
		# Get the deadtime fraction.
		deadTimeHist = rootFile.Get("deadTimeHistograms/hScalarDeadTimeFraction_on")
		binnedLiveTime = {} # A dict is dumb, but used in debugging.
		for i in range(deadTimeHist.GetNbinsX()):
			if i == 0:
				continue
			thisBinCenter = deadTimeHist.GetBinCenter(i)
			thisInterval = (thisBinCenter - 5, thisBinCenter + 5)
			binnedLiveTime[thisInterval] = 1 - deadTimeHist.GetBinContent(i)
		
		binnedLiveTimeArr = np.array(list(binnedLiveTime.values()))
		binnedLiveTimeArr = binnedLiveTimeArr[binnedLiveTimeArr !=0] # Drop zero bins.
		averageLiveTime = np.average(binnedLiveTimeArr) # Get average.
		
		
		elevation = np.average(elevations)
		# Zero the times and get their differences.
		timesZ = times - times[0]
		time_diffs = np.diff(timesZ)
		
		# Too lazy to use np.digitize
		f, ax = plt.subplots()
		n, bins, patches = ax.hist(time_diffs, bins=1000, density=True)
		plt.close()
		
		# Center points
		halfWidth = (bins[-2] - bins[-1])/2
		bins = bins[0:-1] - halfWidth
		p0 = [300]
		xData = bins
		yData = n
	
		# Fit
		try:
			popt, pcov = cf(timeIntervalFunction, xData, yData, p0=p0)
		except Exception as e:
			drop = True
			print(e)

		if not all(np.isfinite(popt)): drop = True
		if not all(np.isfinite(pcov)): drop = True

		rate = popt[0]/averageLiveTime
		rateErr = np.sqrt(np.diagonal(pcov))[0]

		return rate, rateErr, times, elevation, drop, binnedLiveTime



class mscwFileHolder:

	''' For opening many mscw files and getting their event rates. '''
	
	
	def __init__(self, runlistMember, pool, verbosity = 0):
		
		
		self.pool = pool
		self.verbosity = verbosity
		
		self.runlistMember = runlistMember
		
		with open(self.runlistMember.general['runlistFile'], 'r') as f:
			self.runs = [x.strip('\n') for x in f]
			
		columns = [
			'Rate',
			'RateErr',
			'Times', 
			'Elevation',
			'Drop',
			'Zenith',
			'AboveLowerBound',
			'Source',
			'Weather',
			'WeatherNum',
			'ObservationDate',
			'Clean',
			'TimeMask',
			'BinnedData',
			'BinnedDataErr',
			'ObservingSeason',
			'ObservingEpoch',
			'BinnedLiveTime',
			'FitDistance',
			'FitVector',
			'FitParameters'
			]
			
		self.runsCollected = pd.DataFrame(columns=columns)

		self.runsSQLCall = '(' + ','.join(self.runs) + ')'
		self.sourceCall = '''mysql -u readonly -h romulus.ucsc.edu -D VERITAS -r -s -e "select source_id from tblRun_Info where run_id in %s;"''' % self.runsSQLCall 
		self.weatherCall = '''mysql -u readonly -h romulus.ucsc.edu -D VERITAS -r -s -e "select weather from tblRun_Info where run_id in %s;"''' % self.runsSQLCall 
		self.timeCall = '''mysql -u readonly -h romulus.ucsc.edu -D VERITAS -r -s -e "select db_start_time from tblRun_Info where run_id in %s;"''' % self.runsSQLCall

		self.sources = None
		self.weather = None
				
		self.cleanIndx = None
		self.runsClean = None
		self.ratesClean = None
		self.ratesErrClean = None
		self.elevationsClean = None
		self.weatherClean = None
		self.poptFit = None
		self.pcovFit = None
		self.previousClean = None
		self.fitIterations = 0
		
		self.xFill = None
		self.lowerBound = None
		self.upperBound = None
		self.xFit = None
		self.yFit = None
		
		# WHY NOT ADD SOME MORE SHIT
		self.xFillSeasons = {}
		self.lowerBoundSeasons = {}
		self.upperBoundSeasons = {}
		self.xFitSeasons = {}
		self.yFitSeasons = {}	
		self.fitIterationsSeasons = {}
		self.fitParamsSeason = {}
				
	def extractEventRates(self):
		
		sourcesRV = subprocess.Popen(self.sourceCall, stdout=subprocess.PIPE, shell=True).communicate()[0].decode('ascii').split('\n')[:-1]
		weatherRV = subprocess.Popen(self.weatherCall, stdout=subprocess.PIPE, shell=True).communicate()[0].decode('ascii').split('\n')[:-1]
		timeRV = subprocess.Popen(self.timeCall, stdout=subprocess.PIPE, shell=True).communicate()[0].decode('ascii').split('\n')[:-1]
		weatherNums = np.array(self.weatherGrade(weatherRV))
		timeInput = [x.replace(' ', 'T') for x in timeRV]
		observations = Time(timeInput)

		for i, run in enumerate(self.runs):
			
			filename = self.runlistMember.general['saveDir'] + run + '.mscw.root'
			
			thisMSCW = mscwFile(filename, self.verbosity)
			thisMSCW.load()
			if thisMSCW.fileOpen is False: continue
			thisMSCW.extractEventRate()
			thisRate, thisRateErr, theseTimes, thisElevation, thisDrop, thisBinnedLiveTime = thisMSCW.getRateInfo()
			
			runEntry = {
				'Rate' : thisRate,
				'RateErr' : thisRateErr,
				'Times' : theseTimes,
				'Elevation' : thisElevation,
				'Drop' : thisDrop,
				'Zenith' : 90 - thisElevation,
				'AboveLowerBound' : True,
				'Source' : sourcesRV[i],
				'Weather' : weatherRV[i],
				'WeatherNum' : weatherNums[i],
				'ObservationDate' : observations[i],
				'Clean' : False,
				'TimeMask' : [],
				'BinnedData' : [],
				'BinnedDataErr' : [],
				'ObservingSeason' : '',
				'ObservingEpoch' : '',
				'BinnedLiveTime' : thisBinnedLiveTime,
				'FitDistance' : 0.,
				'FitVector' : [0, 0],
				'FitParameters' : []
				}
			

			try:
				self.runsCollected.loc[run] = pd.Series(runEntry)
			except:
				return (runEntry, pd.Series(runEntry))

			if self.verbosity > 1: print('Run ', filename, ' has duration ', (theseTimes.max() - theseTimes.min())/60, ' minutes.')
		
		#print(self.runsCollected)
		seasonsDF = generateSeasonDF()
		for indx, row in seasonsDF.iterrows():
			#print(self.runsCollected['ObservationDate'])
			#print(row['Dates'][0])
			fillIndx = (self.runsCollected['ObservationDate'] >= row['Dates'][0]) & (self.runsCollected['ObservationDate'] <= row['Dates'][1])
			season = row.name
			self.runsCollected.loc[fillIndx, 'ObservingSeason'] = season
		
		#self.previousClean = self.runsCollected['Clean'].copy()
			
	def resetFittedValues(self):
		x = len(self.runsCollected.index)
		self.runsCollected = self.runsCollected.assign(Clean=False)
		self.runsCollected = self.runsCollected.assign(TimeMask=[[] for i in range(x)])
		self.runsCollected = self.runsCollected.assign(BinnedData=[[] for i in range(x)])
		self.runsCollected = self.runsCollected.assign(BinnedDataErr=[[] for i in range(x)])
		self.runsCollected = self.runsCollected.assign(ObservingSeason='')
		self.runsCollected = self.runsCollected.assign(ObservingEpoch='')
		self.runsCollected = self.runsCollected.assign(FitDistance=0.)
		self.runsCollected = self.runsCollected.assign(FitVector=[[0., 0.] for i in range(x)])
		self.runsCollected = self.runsCollected.assign(FitVector=[[] for i in range(x)])
		
	def reseteExtractedValues(self):
		x = len(self.runsCollected.index)
		self.runsCollected = self.runsCollected.assign(Rate=0.)
		self.runsCollected = self.runsCollected.assign(RateErr=0.)
		self.runsCollected = self.runsCollected.assign(Times=[[] for i in range(x)])
		self.runsCollected = self.runsCollected.assign(Elevation=0.)
		self.runsCollected = self.runsCollected.assign(Drop=False)
		self.runsCollected = self.runsCollected.assign(Zenith=0.)
		self.runsCollected = self.runsCollected.assign(AboveLowerBound=True)
		self.runsCollected = self.runsCollected.assign(Source='')
		self.runsCollected = self.runsCollected.assign(Weather='')
		self.runsCollected = self.runsCollected.assign(WeatherNum=0.)	
		self.runsCollected = self.runsCollected.assign(ObservationDate=Time('2001-01-01T00:00:1.000'))	
		self.runsCollected = self.runsCollected.assign(BinnedLiveTime=[{} for i in range(x)])	
	
		
	def resetAllComputedValues(self):
		self.resetFittedValues()
		self.reseteExtractedValues()
		self.sources = None
		self.weather = None
		self.cleanIndx = None
		self.runsClean = None
		self.ratesClean = None
		self.ratesErrClean = None
		self.elevationsClean = None
		self.weatherClean = None
		self.poptFit = None
		self.pcovFit = None
		self.previousClean = None
		self.fitIterations = 0
		self.xFill = None
		self.lowerBound = None
		self.upperBound = None
		self.xFit = None
		self.yFit = None
		self.xFillSeasons = {}
		self.lowerBoundSeasons = {}
		self.upperBoundSeasons = {}
		self.xFitSeasons = {}
		self.yFitSeasons = {}	
		self.fitIterationsSeasons = {}
		self.fitParamsSeason = {}

		
	@staticmethod
	def weatherGrade(skys):
		
		skyOut = []
		
		for sky in skys:
			skyNum = -1
			if sky[0] == 'A':
				skyNum = 4
			elif sky[0] == 'B':
				skyNum = 3
			elif sky[0] == 'C':
				skyNum = 2
			elif sky[0] == 'D':
				skyNum = 1
			elif sky[0] == 'F':
				skyNum = 0
			try:
				if sky[1] == '+':
					skyNum += 0.3
				elif sky[1] == '-':
					skyNum -= 0.3
			except:
				pass
			
			skyOut.append(skyNum)
			
		return skyOut
	
	def binRuns(self):
		for indx, row in self.runsCollected.iterrows():
			binningInterval = self.runlistMember.parameters['binningInterval']
			binnedRate, binnedRateErr = binRate(row['Times'], self.pool, interval=binningInterval)
			nLivetimes = len(list(row['BinnedLiveTime'].values()))
			nBinnedRates = len(binnedRate)
			#print('nLivetimes', nLivetimes)
			#print('nBinnedRates', nBinnedRates)
			if nLivetimes == nBinnedRates:
				binnedRate = binnedRate/np.array(list(row['BinnedLiveTime'].values()))
			elif nBinnedRates > nLivetimes:
				nDiff = nBinnedRates - nLivetimes
				thisLivetime = np.pad(np.array(list(row['BinnedLiveTime'].values())), nDiff, 'constant', constant_values=1)
				binnedRate = binnedRate/thisLivetime
			elif nBinnedRates < nLivetimes:
				binnedRate = binnedRate/np.array(list(row['BinnedLiveTime'].values()))[0:nBinnedRates]
			
			self.runsCollected.loc[indx, 'BinnedData'].append(binnedRate)
			self.runsCollected.loc[indx, 'BinnedDataErr'].append(binnedRateErr)

	def plotRun(self, run):
		try:
			indx = np.where(self.runs == str(run))[0]
		except Exception as e:
			print('Failed to locate run. Returning.')
			print(e)
			return
		
		filename = self.runlistMember.general['saveDir'] + str(run) + '.mscw.root'
		
		thisMSCW = mscwFile(filename, self.verbosity)
		thisMSCW.load()
		thisMSCW.extractEventRate()
		thisRate, thisRateErr, theseTimes, thisElevation, thisDrop = thisMSCW.getRateInfo()
		print('Run Values:')
		print('Run: ', run)
		print('Rate: ', thisRate)
		print('Rate Error: ', thisRateErr)
		print('thisElevation: ', thisElevation)
		time_diffs = np.diff(theseTimes)
	
		f, ax = plt.subplots()
		ax.set_xlabel('Time difference [s]')
		ax.set_ylabel('I(t)')
		ax.set_title(str(run))
		
		n1, bins1, patches1 = ax.hist(time_diffs, bins=1000, density=True, label=str(run) + ' Time interval distribution', alpha=0.75)
		halfWidth1 = (bins1[-2] - bins1[-1])/2
		bins1 = bins1[0:-1] - halfWidth1
	
	
		p0 = [300]
		xData1 = bins1
		yData1 = n1
	
		try:
			popt1, pcov1 = cf(timeIntervalFunction, xData1, yData1, p0=p0)
		except Exception as e:
			print(e)
			print('Failed to get event rate on run ', run)
		
		xFit1 = np.linspace(0, np.max(xData1), 1000)
		yFit1 = timeIntervalFunction(xFit1, *popt1)
		ax.plot(xFit1, yFit1)
		
		plt.show(block=True)
		

		
		
	def plotTimeMasks(self, show=True):
		for indx, row in self.runsCollected.iterrows():
			
			#print('Working on run ', indx)
			stdFactor = self.runlistMember.parameters['stdFactor']
			
			binnedRate = row['BinnedData'][0]
			binnedRateErr = row['BinnedDataErr'][0]
			medianRate = np.average(binnedRate)
			stdRate = np.std(binnedRate)
			startTime = row['Times'][0]
			endTime = row['Times'][-1]
			duration = endTime - startTime # Seconds
			roundSeconds = int(duration)	

			zoomFactor = float(roundSeconds/binnedRate.shape[0])
			zoomedRates = ni.zoom(binnedRate, zoomFactor)
			zeroedRates = zoomedRates - medianRate

			fig, ax = plt.subplots()
			ax.plot(np.linspace(0, duration, binnedRate.shape[0]), binnedRate, color='black')
			ax.plot(np.linspace(0, duration, zoomedRates.shape[0]), zoomedRates, color='red')
			#ax.scatter(np.linspace(0, duration, row['BinnedLiveTime'].shape[0]), row['BinnedLiveTime'], color='blue')
			#ax.plot(np.linspace(0, duration, 2), [-1 * stdFactor * stdRate, -1 * stdFactor * stdRate])
			ax.set_xlabel('Time (30 s intervals) [s]')
			ax.set_ylabel('Binned rate [Hz]')
			
			for mask in row['TimeMask']:
				ax.fill_between(mask, [-100, -100], [100, 100], color='black', alpha=0.25)
				
			#ax.set_xlim([0, duration])
			#ax.set_ylim([zeroedRates.min(), zeroedRates.max()])
			extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
			runInfo = 'Run ' + indx
			stdInfo = 'Std ' + str(stdRate)
			rateInfo = 'Rate ' + str(medianRate)
			ax.legend([extra, extra, extra], [runInfo, stdInfo, rateInfo])
			ax.set_title(indx)

			if self.runlistMember.general['savePlots'] is True:
				subprocess.call('mkdir -p ' + self.runlistMember.general['diagnosticFolder'],shell=True)
				plt.savefig(self.runlistMember.general['diagnosticFolder'] + 'Run_' + indx + '_timemask_cut.png')
			if show == True:
				plt.show()

			plt.close()
			

	def calcTimeMaskParameters(self, row):
		
		binnedRate = row['BinnedData'][0]
		binnedRateErr = row['BinnedDataErr'][0]
		medianRate = np.median(binnedRate)
		avgRate = np.average(binnedRate)
		stdRate = np.std(binnedRate)
		startTime = row['Times'][0]
		endTime = row['Times'][-1]
		duration = endTime - startTime # Seconds
		roundSeconds = int(duration)
		timeInfo = (startTime, endTime, duration, roundSeconds)
		
		return (binnedRate, binnedRateErr, medianRate, avgRate, stdRate, timeInfo)
	
	def calcIntrepRate(self, binnedRate, roundSeconds):
		
		zoomFactor = float(roundSeconds/binnedRate.shape[0])
		zoomedRates = ni.zoom(binnedRate, zoomFactor)
		return zoomedRates
			
	def resetTimeMasks(self):
		
		for indx, row in self.runsCollected.iterrows():
			row['TimeMask'] = []
			
	def plotTimeMasksFunc(self, indx, show=False, save=False):
		
		medianRate = self.runsCollected.loc[indx, 'medianRate']
		stdRate = self.runsCollected.loc[indx, 'stdRate']
		stdFactor = self.runsCollected.loc[indx, 'stdFactor']
		duration = self.runsCollected.loc[indx, 'duration']
		zeroedRate = self.runsCollected.loc[indx, 'zeroedRate']			
		
		fig, ax = plt.subplots()
		ax.set_xlabel('Time (30 s intervals) [s]')
		ax.set_ylabel('Binned rate [Hz]')
		
		
		stdThresh = -1 * stdFactor * stdRate
		#ax.plot(np.linspace(0, duration, binnedRate.shape[0]), binnedRate, color='black')
		ax.plot(np.linspace(0, duration, zeroedRate.shape[0]), zeroedRate, color='black')
		ax.plot([-10, duration + 10], [stdThresh, stdThresh], color='red', ls='--')

		for mask in self.runsCollected.loc[indx, 'TimeMask']:
			ax.fill_between(mask, [-100, -100], [100, 100], color='black', alpha=0.25)
			
		ax.set_xlim([0, duration])
		ax.set_ylim([zeroedRate.min(), zeroedRate.max()])
		extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
		runInfo = 'Run ' + indx
		stdInfo = 'Std ' + str(stdRate)
		rateInfo = 'Rate ' + str(medianRate)
		ax.legend([extra, extra, extra], [runInfo, stdInfo, rateInfo])
		ax.set_title(indx)

		if save is True:
			subprocess.call('mkdir -p ' + self.runlistMember.general['diagnosticFolder'],shell=True)
			plt.savefig(self.runlistMember.general['diagnosticFolder'] + 'Run_' + indx + '_timemask_cut.png')
		if show is True:
			plt.show()

		plt.close()		
			
	def intraRunMask(self, plot=True, show=False, median=False):
				
		minCutIntervals = []
		stdFactors = []
		paddingTimes = []
		combineSeperations = []
		endThresholds = []
		binnedRates = []
		binnedRateErrs = []
		medianRates = []
		avgRates = []
		stdRates = []
		startTimes = []
		endTimes = []
		durations = []
		roundSecondss = []
		zoomedRatess = []
		zeroedRatess = []
		
		for indx, row in self.runsCollected.iterrows():
			#print('Working on run ', indx)
			
			minCutInterval = self.runlistMember.parameters['minCutInterval']
			stdFactor = self.runlistMember.parameters['stdFactor']
			paddingTime = self.runlistMember.parameters['paddingTime']
			combineSeperation = self.runlistMember.parameters['combineSeperation']
			endThreshold = self.runlistMember.parameters['endThreshold']

			binnedRate, binnedRateErr, medianRate, avgRate, stdRate, timeInfo = self.calcTimeMaskParameters(row)
			startTime, endTime, duration, roundSeconds = timeInfo
			zoomedRates = self.calcIntrepRate(binnedRate, roundSeconds) # Lazy; go from 10 secs bins -> 1 sec by interp.

			
			if median == False:
				zeroedRates = zoomedRates - avgRate
			else:
				zeroedRates = zoomedRates - medianRate
				
			minCutIntervals.append(minCutInterval)
			stdFactors.append(stdFactor)
			paddingTimes.append(paddingTime)
			combineSeperations.append(combineSeperation)
			endThresholds.append(endThreshold)
			binnedRates.append(binnedRate)
			binnedRateErrs.append(binnedRateErr)
			medianRates.append(medianRate)
			avgRates.append(avgRate)
			stdRates.append(stdRate)
			startTimes.append(startTime)
			endTimes.append(endTime)
			durations.append(duration)
			roundSecondss.append(roundSeconds)
			zoomedRatess.append(zoomedRates)
			zeroedRatess.append(zeroedRates)
			
			if [0, 10000] in row['TimeMask']:
				if self.verbosity > 2: print('Skipping run ', indx, ', it is marked to wholly drop.')
				continue
			
			if stdRate < self.runlistMember.parameters['stdThresh']:
				if self.verbosity > 2:
					print('Run ', indx, ' standard deviation below threshold, not doing intrarun timemasking. ')
				continue
			
			negativeRates = np.where(zeroedRates <= -1 * stdFactor * stdRate)[0]
			regions = consecutive(negativeRates)
			regions = [x for x in regions if len(x) > minCutInterval]
			if len(regions) < 1:
				if self.verbosity > 2:
					print('Run ', indx, ' had no regions above minimum cut interval and standard deviation threshold, no intrarun timemask. ')
				continue
			timeIntervals = [[max(region[0]-paddingTime, 0), region[-1]+paddingTime] for region in regions]
			
			intervals = len(timeIntervals)
			combinedIntervals = []
			if len(timeIntervals) > 1:
				for j, timeInterval in enumerate(timeIntervals):
					if j < intervals - 1:
						intervalStart = timeIntervals[j][-1]
						intervalEnd = timeIntervals[j+1][0]
						intervalGap = intervalEnd - intervalStart
						if intervalGap > combineSeperation:
							combinedIntervals.append(timeIntervals[j])
						else:
							combinedInterval = [timeIntervals[j][0], timeIntervals[j+1][-1]]
							combinedIntervals.append(combinedInterval)
					else:
						intervalEnd = timeIntervals[j][0]
						intervalStart = timeIntervals[j-1][-1]
						intervalGap = intervalEnd - intervalStart
						if intervalGap > combineSeperation:
							combinedIntervals.append(timeIntervals[j])
			else:
				combinedIntervals = timeIntervals
				
			if combinedIntervals[-1][-1] + endThreshold > duration:
				combinedIntervals[-1][-1] = duration + 1

			totalMaskTime = np.array([x[-1] - x[0] for x in combinedIntervals]).sum()
				
			if totalMaskTime > 0.75 * duration or stdRate > 20:
				if self.verbosity > 2:
					print('Run ', indx, ' had total time mask fraction above threshold or standard deviation too high, masking whole run.')
				combinedIntervals = [[0, 10000]]
				
			for ci in combinedIntervals:
				self.runsCollected.loc[indx, 'TimeMask'].append(ci)
				
			if self.verbosity > 1:
				print('++++++++++++++++++++++++++++++')
				print('Run: ', indx)
				print('Duration: ', duration)
				print('std: ', stdRate)
				print('Intrarun mask(s): ', combinedIntervals)
				print('++++++++++++++++++++++++++++++')
				
			print(minCutIntervals)
			
		#import pdb; pdb.set_trace()
			
		self.runsCollected['MinCutInterval'] = minCutIntervals
		self.runsCollected['stdFactor'] = stdFactors
		self.runsCollected['paddingTime'] = paddingTimes
		self.runsCollected['combineSeperation'] = combineSeperations
		self.runsCollected['endThreshold'] = endThresholds
		self.runsCollected['medianRate'] = medianRates
		self.runsCollected['avgRate'] = avgRates
		self.runsCollected['stdRate'] = stdRates
		self.runsCollected['startTime'] = startTimes
		self.runsCollected['endTime'] = endTimes
		self.runsCollected['duration'] = durations
		self.runsCollected['roundSeconds'] = roundSecondss
		self.runsCollected['zoomedRate'] = zoomedRatess
		self.runsCollected['zeroedRate'] = zeroedRatess
		
		if plot is True:
			for indx, row in self.runsCollected.iterrows():
				self.plotTimeMasksFunc(indx, show=show, save=self.runlistMember.general['savePlots'])

	def intraRunMask2(self, plot=True, show=False):
				
		minCutIntervals = []
		stdFactors = []
		paddingTimes = []
		combineSeperations = []
		endThresholds = []
		binnedRates = []
		binnedRateErrs = []
		medianRates = []
		avgRates = []
		stdRates = []
		startTimes = []
		endTimes = []
		durations = []
		roundSecondss = []
		zoomedRatess = []
		zeroedRatess = []
		
		for indx, row in self.runsCollected.iterrows():

			minCutInterval = self.runlistMember.parameters['minCutInterval']
			stdFactor = self.runlistMember.parameters['stdFactor']
			paddingTime = self.runlistMember.parameters['paddingTime']
			combineSeperation = self.runlistMember.parameters['combineSeperation']
			endThreshold = self.runlistMember.parameters['endThreshold']
			

			
			#if indx != '67227':
			#	continue
			duration = row['Times'][-1] - row['Times'][0]
			binnedRate, binnedRateErr = binRate(row['Times'], self.pool, interval=1)
			zoomedLivetime = getInterpLiveTime(np.array(list(row['BinnedLiveTime'].values())), binnedRate.shape[0])
			binnedRate = binnedRate/zoomedLivetime
			binnedRateErr = binnedRateErr/zoomedLivetime
			window = signal.general_gaussian(20, p=1, sig=5)
			smoothedData = signal.fftconvolve(window, binnedRate, mode='full')
			smoothedData = (np.average(binnedRate) / np.average(smoothedData)) * smoothedData
			smoothedData = np.roll(smoothedData, -1*int((20-1)/2))
			smoothedData[0:int(20/2)] = smoothedData[int(20/2)]
			smoothedData[-int(60/2):] = smoothedData[-int(60/2) -1]
			

			
			medianRate = np.median(smoothedData)
			avgRate = np.average(smoothedData)
			startTime = row['Times'][0]
			endTime = row['Times'][-1]
			roundSeconds = int(duration)
			zoomedRates = smoothedData
			
			
			zeroedData = smoothedData - medianRate
			swf = splineCorrection2(zeroedData, finalArrSize=binnedRate.shape[0], nDivisions=4)
			excess = np.sum(np.abs(swf))/swf.shape[0]
			zeroedRates = zeroedData
			
			stdRate = np.std(zeroedData[np.where(zeroedData > 0)[0]])
			
			stdFactor = 3
		
			stdThresh = -1 * stdFactor * stdRate
			
			minCutIntervals.append(minCutInterval)
			stdFactors.append(stdFactor)
			paddingTimes.append(paddingTime)
			combineSeperations.append(combineSeperation)
			endThresholds.append(endThreshold)
			binnedRates.append(binnedRate)
			binnedRateErrs.append(binnedRateErr)
			medianRates.append(medianRate)
			avgRates.append(avgRate)
			stdRates.append(stdRate)
			startTimes.append(startTime)
			endTimes.append(endTime)
			durations.append(duration)
			roundSecondss.append(roundSeconds)
			zoomedRatess.append(zoomedRates)
			zeroedRatess.append(zeroedRates)
			
			negativeRates = np.where(zeroedData <= stdThresh)[0]
			regions = consecutive(negativeRates)
			regions = [x for x in regions if len(x) > minCutInterval]
			if len(regions) < 1:
				continue
			timeIntervals = [[max(region[0]-paddingTime, 0), region[-1]+paddingTime] for region in regions]
			
			intervals = len(timeIntervals)
			combinedIntervals = []
			if len(timeIntervals) > 1:
				for j, timeInterval in enumerate(timeIntervals):
					if j < intervals - 1:
						intervalStart = timeIntervals[j][-1]
						intervalEnd = timeIntervals[j+1][0]
						intervalGap = intervalEnd - intervalStart
						if intervalGap > combineSeperation:
							combinedIntervals.append(timeIntervals[j])
						else:
							combinedInterval = [timeIntervals[j][0], timeIntervals[j+1][-1]]
							combinedIntervals.append(combinedInterval)
					else:
						intervalEnd = timeIntervals[j][0]
						intervalStart = timeIntervals[j-1][-1]
						intervalGap = intervalEnd - intervalStart
						if intervalGap > combineSeperation:
							combinedIntervals.append(timeIntervals[j])
			else:
				combinedIntervals = timeIntervals
				
						
			if combinedIntervals[-1][-1] + endThreshold > duration:
				combinedIntervals[-1][-1] = duration + 1

			totalMaskTime = np.array([x[-1] - x[0] for x in combinedIntervals]).sum()
				
			if totalMaskTime > 0.75 * duration or stdRate > 20:
				if self.verbosity > 1:
					print('Run ', indx, ' had total time mask fraction above threshold or standard deviation too high, masking whole run.')
				combinedIntervals = [[0, 10000]]

			if excess > 20:
				if self.verbosity > 1:
					print('Run ', indx, ' had spline excess above threshold, masking whole run.')
				combinedIntervals = [[0, 10000]]	

			if [0, 10000] in row['TimeMask']:
				if self.verbosity > 1: print('Skipping run ', indx, ', it is marked to wholly drop.')
				combinedIntervals = [[0, 10000]]
				
			for ci in combinedIntervals:
				self.runsCollected.loc[indx, 'TimeMask'].append(ci)
				
			if self.verbosity > 2:
				print('++++++++++++++++++++++++++++++')
				print('Run: ', indx)
				print('Duration: ', duration)
				print('std: ', stdRate)
				print('Intrarun mask(s): ', combinedIntervals)
				print('++++++++++++++++++++++++++++++')

			f, ax = plt.subplots()

			ax.plot(zeroedData, color='blue', label='Smoothed L3 rate')
			ax.plot([-10, duration + 10], [stdThresh, stdThresh], color='red', ls='--', label='Threshold')
			ax.scatter(negativeRates, zeroedData[negativeRates], marker='.', color='black', label='Points below threshold')
			ax.plot(swf, color='Green', label='Spline fit, excess=%2.2f' % (excess))
			ax.set_title('Run ' + indx)
			ax.set_xlabel('Time [s]')
			ax.set_ylabel('Frequency (zeroed) [Hz]')
			ax.set_xlim([0, duration])
			ax.set_ylim([zeroedData.min()-20, zeroedData.max()+20])
			
			for ci in combinedIntervals:
				ax.fill_between(ci, [-1000, -1000], [1000, 1000], color='black', alpha=0.25)
			plt.legend()	

				
			if plot is True:
				plt.savefig(self.runlistMember.general['diagnosticFolder'] + 'Run_' + indx + '_timemask_cut.png')
			if show is True:
				plt.show()
			plt.close()
				
				
			#print(minCutIntervals)
			
		#import pdb; pdb.set_trace()
			
		self.runsCollected['MinCutInterval'] = minCutIntervals
		self.runsCollected['stdFactor'] = stdFactors
		self.runsCollected['paddingTime'] = paddingTimes
		self.runsCollected['combineSeperation'] = combineSeperations
		self.runsCollected['endThreshold'] = endThresholds
		self.runsCollected['medianRate'] = medianRates
		self.runsCollected['avgRate'] = avgRates
		self.runsCollected['stdRate'] = stdRates
		self.runsCollected['startTime'] = startTimes
		self.runsCollected['endTime'] = endTimes
		self.runsCollected['duration'] = durations
		self.runsCollected['roundSeconds'] = roundSecondss
		self.runsCollected['zoomedRate'] = zoomedRatess
		self.runsCollected['zeroedRate'] = zeroedRatess
		
		#if plot is True:
		#	for indx, row in self.runsCollected.iterrows():
		#		self.plotTimeMasksFunc(indx, show=show, save=bool(self.runlistMember.general['savePlots']))


	def callSQLDatabase(self):
		
		self.sourcesRV = subprocess.Popen(self.sourceCall, stdout=subprocess.PIPE, shell=True).communicate()[0].decode('ascii').split('\n')[:-1]
		self.weatherRV = subprocess.Popen(self.weatherCall, stdout=subprocess.PIPE, shell=True).communicate()[0].decode('ascii').split('\n')[:-1] 
		
		self.sources = set(self.sourcesRV)
		self.weather = np.array(self.weatherGrade(self.weatherRV))
		
	def plotRuns(self, markOutliers = False, markAll = False, show=True):

		cleanRuns = self.runsCollected['Clean']		
		cleanZeniths = self.runsCollected['Zenith'].loc[cleanRuns]
		cleanRates = self.runsCollected['Rate'].loc[cleanRuns]
		
		if self.verbosity > 5:
			print('Using this dataframe:')
			print(self.runsCollected)
		
		f, ax = plt.subplots(figsize=(16,8))
		#ax.set_aspect('equal')
		rv = ax.scatter(cleanZeniths, cleanRates, marker='.', s=750, label=r'Weather-rated runs', color='yellow', cmap=cm.plasma, edgecolors='black', alpha=0.5)
		rv = ax.scatter(self.runsCollected['Zenith'], self.runsCollected['Rate'], marker='.', s=500, label=r'Weather-rated runs', c=self.runsCollected['WeatherNum'], cmap=cm.plasma, edgecolors='black')
		cb = f.colorbar(rv)
		cb.set_label('Sky Rating')
		ax.plot(self.xFit, self.yFit, label='Cosine Fit')
		plt.fill_between(self.xfill, self.lowerBound, self.upperBound*100, color = 'black', alpha = 0.15, label='10% Band')
		max1 = np.max(self.runsCollected['Rate']) + 20
		max2 = np.max(self.yFit + 20)
		maxs = np.array([max1, max2])
		ax.set_ylim(0, np.max(maxs))
		ax.set_xlim(0, 90)
		ax.set_xlabel(r'Zenith Angle [$^\circ$]')
		ax.set_ylabel('Event rate after size cut [Hz]')
		#ax.set_title(self.runlistMember.parameters['saveNameStem'] + ' popt = ' + str(self.poptFit))
		countObject = str(Counter(self.runsCollected['Source']).most_common(1)[0][0])
		countSeason = str(Counter(self.runsCollected['ObservingSeason']).most_common(1)[0][0])
		#ax.set_title(countObject +  ' ' + countSeason + ' season')
		filename = countObject +  '_' + countSeason + '_season.png'
		#plt.legend()
		if markAll is True:
			for indx, row in self.runsCollected.iterrows():
				ax.text(row['Zenith'], row['Rate'], indx + '\n' + ' (%f, %f) ' % (round(row['Zenith'], 2), round(row['Rate'], 2)), fontsize = 8)		
		if markOutliers is True:
			for indx, row in self.runsCollected.iterrows():
				if row['AboveLowerBound'] is False:
					ax.text(row['Zenith'], row['Rate'], indx + '\n' + ' (%f, %f) ' % (round(row['Zenith'], 2), round(row['Rate'], 2)), fontsize = 8)		
		#print('Popt:')
		#print(self.poptFit)
		#print('Pcov:')
		#print(self.pcovFit)
		
		plt.tight_layout()
		if self.runlistMember.general['savePlots'] is True:
			subprocess.call('mkdir -p ' + self.runlistMember.general['diagnosticFolder'],shell=True)
			plt.savefig(self.runlistMember.general['diagnosticFolder'] + filename)
		if show is True:
			plt.show()
		plt.close()

		
	def excludeRuns(self):
		for indx, row in self.runsCollected.iterrows():
			if row['AboveLowerBound'] is False:
				self.runsCollected.loc[indx, 'TimeMask'].append([0, 10000])
			

					
	def writeTimeMask(self):
		with open(os.path.expandvars(self.runlistMember.parameters['timemaskFile']), 'w') as f:
			for indx, row in self.runsCollected.iterrows():
				for mask in row['TimeMask']:
					f.write('* ' + indx + ' %d %d 0 \n' % (mask[0], mask[1]))
		
	def identifyOutliers(self):
		
		if self.lowerBound is not None and self.xfill is not None:
			for indx, row in self.runsCollected.iterrows():
				
				thisX = row['Zenith']
				thisY = row['Rate']
				checkX = self.xfill[np.argmin(np.abs(self.lowerBound - thisY))]
				checkY = self.lowerBound[np.argmin(np.abs(self.xfill - thisX))]
				distPoint = np.sqrt(thisX**2 + thisY**2)
				distCheck = np.sqrt(checkX**2 + checkY**2)
				
				
				#if thisX < checkX and thisY < checkY:
				if distPoint < distCheck:
					self.runsCollected.loc[indx, 'AboveLowerBound'] = False
		else:
			print('Lower bound threshold not set, run the fit.')

		
	
	def findCleanRuns(self):
		
		self.runsCollected['Clean'] = False
		self.runsCollected.loc[(self.runsCollected['WeatherNum'] >= 4.0) & (self.runsCollected['Zenith'] <= 80) & (self.runsCollected['AboveLowerBound']), 'Clean'] = True
		#cleanRuns = (self.runsCollected['Clean'] == True)
		#cleanZeniths = self.runsCollected.loc[cleanRuns, 'Zenith']
		#if len(cleanZeniths) == 0:
		#	self.runsCollected.loc[(self.runsCollected['WeatherNum'] >= 3.0) & (self.runsCollected['Zenith'] <= 80) & (self.runsCollected['AboveLowerBound']), 'Clean'] = True

			
	def iterateFit(self, plot=False, **kwargs):
		
		
		seasons = set(self.runsCollected['ObservingSeason'])
		self.runsCollectedOrig = self.runsCollected.copy()
		self.runsCollectedTemp = self.runsCollected.copy()
		
		for s, season in enumerate(seasons):
			
			self.fitIterations = 0
			thisSeason = self.runsCollectedOrig['ObservingSeason'] == season
			self.runsCollected = self.runsCollectedOrig.loc[thisSeason].copy()
			self.previousClean = self.runsCollected['Clean'].copy()
			
			#print('On this season: ', season)
			#print('With this df: ', self.runsCollected)
			
			if np.where(self.runsCollected['WeatherNum'] >= 4.0)[0].shape[0] == 0:
				print('Not enough clean weather runs to fit this season. Skipping.')
				continue

			while not ((self.runsCollected['Clean'] == self.previousClean).all() and self.fitIterations != 0):
				self.previousClean = self.runsCollected['Clean'].copy()
				self.findCleanRuns()
				self.fit()
				self.identifyOutliers()
				self.fitIterations += 1
				self.calculateDistances2()
				
				if plot is True:
					self.plotRuns(**kwargs)
				if self.fitIterations > 10:
					print('Iterations greater than 10, returning.')
					continue
				
			self.xFillSeasons[season] = self.xFill
			self.lowerBoundSeasons[season] = self.lowerBound
			self.upperBoundSeasons[season] = self.upperBound
			self.xFitSeasons[season] = self.xFit
			self.yFitSeasons[season] = self.yFit	
			self.fitIterationsSeasons[season] = self.fitIterations
			
			self.runsCollectedTemp[thisSeason] = self.runsCollected
		self.runsCollectedOrig = self.runsCollectedTemp
		self.runsCollected = self.runsCollectedTemp
		
		#bounds = ((50, 0), (500, 10))
		# ,   0.91702798 p0
		# 0 10 bounds
	def fit(self, p0 = [146.96540466], bounds = ((50), (500)), thres = 0.10):

		cleanRuns = (self.runsCollected['Clean'] == True)
		cleanZeniths = self.runsCollected.loc[cleanRuns, 'Zenith']
		cleanRates = self.runsCollected.loc[cleanRuns, 'Rate']
		cleanRateErr = self.runsCollected.loc[cleanRuns, 'RateErr']
		
		if len(cleanZeniths) == 0:
			print('WARNING! NO CLEAN RUNS TO FIT. BAD VALUE!')
			self.poptFit = np.array([50.])
			self.xFit = np.linspace(0, 90, 1000)
			self.yFit = cosineFit(self.xFit, *self.poptFit)
			self.pcovFit = np.array([[np.inf]])
			
			xs = self.xFit[1:]
			ys = self.yFit[1:]
	
			u = np.diff(self.xFit)
			v = np.diff(self.yFit)
	
			u20 = (thres*self.yFit[0])*u/np.sqrt(u**2 + v**2)
			v20 = (thres*self.yFit[0])*v/np.sqrt(u**2 + v**2)	
	
			xTipsLower20 = xs + v20
			xTipsUpper20 = xs + -1*v20
			yTipsLower20 = ys + -1*u20
			yTipsUpper20 = ys + u20
	
			self.xfill = np.sort(np.concatenate([xTipsLower20, xTipsUpper20]))
			self.lowerBound = np.interp(self.xfill, xTipsLower20, yTipsLower20)
			self.upperBound = np.interp(self.xfill, xTipsUpper20, yTipsUpper20)
			self.fitIterations = 100
			return


		if self.verbosity > 3:
			print('Fitting on these ',len(cleanRuns), 'zenith values: \n', self.runsCollected.loc[cleanRuns, ('Zenith', 'Rate', 'Weather', 'WeatherNum')])
		
		self.poptFit, self.pcovFit = cf(cosineFit, cleanZeniths, cleanRates, sigma = cleanRateErr, maxfev=10000, bounds=bounds)
		self.xFit = np.linspace(0, 90, 1000)
		self.yFit = cosineFit(self.xFit, *self.poptFit)
		yFitErr = np.sqrt(np.diagonal(self.pcovFit))
		if self.verbosity > 3:
			print('popt and pcov:')
			print(self.poptFit)
			print(self.pcovFit)

		xs = self.xFit[1:]
		ys = self.yFit[1:]

		u = np.diff(self.xFit)
		v = np.diff(self.yFit)

		u20 = (thres*self.yFit[0])*u/np.sqrt(u**2 + v**2)
		v20 = (thres*self.yFit[0])*v/np.sqrt(u**2 + v**2)	

		xTipsLower20 = xs + v20
		xTipsUpper20 = xs + -1*v20
		yTipsLower20 = ys + -1*u20
		yTipsUpper20 = ys + u20

		self.xfill = np.sort(np.concatenate([xTipsLower20, xTipsUpper20]))
		self.lowerBound = np.interp(self.xfill, xTipsLower20, yTipsLower20)
		self.upperBound = np.interp(self.xfill, xTipsUpper20, yTipsUpper20)
	
		
	def calculateDistances2(self):
	
		
		for i, row in self.runsCollected.iterrows():
			self.runsCollected.loc[i, 'FitParameters'] = self.poptFit
			
			
			thisX = row['Zenith']
			thisY = row['Rate']
			theDistances = distance(self.xFit, self.poptFit, thisX, thisY)
			minDistance = theDistances.min()
			minArg = np.argmin(theDistances)
			minX = self.xFit[minArg]
			minY = self.yFit[minArg]
			self.runsCollected.loc[i, 'FitDistance'] = minDistance
			self.runsCollected.loc[i, 'FitVector'][0] = minX
			self.runsCollected.loc[i, 'FitVector'][1] = minY
			#print(i, ' Distance: ', minDistance)
			#print(i, ' minX: ', minX)
			#print(i, ' minY: ', minY)			
			
		#print(self.runsCollected['FitDistance'])
			

	def plotHistos(self):
	
		def lineFit(x, m, b):
			return m*x + b
		
		seasons = set(self.runsCollected['ObservingSeason'])
		seasonsNorm = np.array([])
		seasonsNum = np.array([])
		seasonsList = np.array(list(seasons))
		summerDists = np.array([])
		winterDists = np.array([])
		
		for s, season in enumerate(seasons):
			
			source = self.runsCollected['Source'][0]
			f, ax = plt.subplots()
			
			theseRuns = self.runsCollected['ObservingSeason'] == season 
			
			dists = self.runsCollected.loc[theseRuns, 'FitDistance']
			if season[-1] == 'W':
				winterDists = np.append(winterDists, dists)
			else:
				summerDists = np.append(summerDists, dists)
			
			print('Season ', season, ': ', np.average(dists))
			cut = 0.1 * self.runsCollected.loc[theseRuns, 'FitParameters'][0]
			
			n, bins, patches = plt.hist(dists, bins = len(np.arange(0, dists.max(), 1)), color= 'black', label='Distance distribution')
			maxY = n.max()
			
			ax.plot([cut, cut], [0, maxY], color='red', ls='--', label='Cut threshold')
			ax.set_title(source + ' ' + season + ' distances-from-fit distribution')
			ax.set_xlabel('Run distance from fit [Hz]')
			ax.set_ylabel('N')
			
			plt.legend()
			filename = source +  '_' + season + '_season_distro.png'
			subprocess.call('mkdir -p ' + self.runlistMember.general['diagnosticFolder'],shell=True)
			plt.savefig(self.runlistMember.general['diagnosticFolder'] + filename)
			plt.close()
			
			seasonsNorm = np.append(seasonsNorm, self.runsCollected.loc[theseRuns, 'FitParameters'][0])
			year = le(season[0:-1])
			mid = 0
			if season[-1] == 'W':
				mid = 0
			else:
				mid = 0.5
			seasonsNum = np.append(seasonsNum, year + mid)
	
	
		f, ax = plt.subplots()
		dists = self.runsCollected['FitDistance']
		
		n, bins, patches = plt.hist(dists, bins = len(np.arange(0, dists.max(), 1)), color= 'black', label='Distance distribution')
		maxY = n.max()
		
		#ax.plot([cut, cut], [0, maxY], color='red', ls='--', label='Cut threshold')
		ax.set_title('All runs, distances-from-fit distribution')
		ax.set_xlabel('Run distance from fit [Hz]')
		ax.set_ylabel('N')
		
		plt.legend()
		filename = source +  '_total_distro.png'
		subprocess.call('mkdir -p ' + self.runlistMember.general['diagnosticFolder'],shell=True)
		plt.savefig(self.runlistMember.general['diagnosticFolder'] + filename)
		plt.close()
		
		##
		
		f, ax = plt.subplots()
		
		n, bins, patches = plt.hist(summerDists, bins = len(np.arange(0, summerDists.max(), 1)), color= 'black', label='Summer Distance Distribution')
		maxY = n.max()
		
		#ax.plot([cut, cut], [0, maxY], color='red', ls='--', label='Cut threshold')
		ax.set_title('Summer runs, distances-from-fit distribution')
		ax.set_xlabel('Run distance from fit [Hz]')
		ax.set_ylabel('N')
		
		plt.legend()
		filename = source +  '_summerTotal_distro.png'
		plt.savefig(self.runlistMember.general['diagnosticFolder'] + filename)
		plt.close()
		
		f, ax = plt.subplots()
		
		n, bins, patches = plt.hist(winterDists, bins = len(np.arange(0, winterDists.max(), 1)), color= 'black', label='Summer Distance Distribution')
		maxY = n.max()
		
		#ax.plot([cut, cut], [0, maxY], color='red', ls='--', label='Cut threshold')
		ax.set_title('All runs, distances-from-fit distribution')
		ax.set_xlabel('Winter distance from fit [Hz]')
		ax.set_ylabel('N')
		
		plt.legend()
		filename = source +  '_winterTotal_distro.png'
		plt.savefig(self.runlistMember.general['diagnosticFolder'] + filename)
		plt.close()
	
		##
		args = np.argsort(seasonsNum)
		seasonsNum = seasonsNum[args]
		seasonsNorm = seasonsNorm[args]
		seasonsList = seasonsList[args]
		
		f, ax = plt.subplots()
	
		ax.scatter(seasonsNum, seasonsNorm, color='black', label='Winter Seasons')
		summers = seasonsNum % 1 == 0.5
		ua = seasonsNum >= 2013
		popt, pcov = cf(lineFit, seasonsNum[ua], seasonsNorm[ua])
		xs = np.linspace(2012, 2025, 10)
		ys = lineFit(xs, *popt)
		ax.plot(xs, ys, color='blue', label='Linear Fit')
		
		xs = np.linspace(2006, 2012, 10)
		ys = lineFit(xs, *popt)
		ax.plot(xs, ys, color='blue', ls='--')
		
		ax.scatter(seasonsNum[summers], seasonsNorm[summers], color='red', label='Summer Seasons')
		ax.set_xlabel('Season')
		ax.set_ylabel('Rate normalization [Hz]')
		ax.set_xlim([2006, 2022.5])
		
		filename = source +  '_scatter.png'
		plt.savefig(self.runlistMember.general['diagnosticFolder'] + filename)
		plt.legend()
		plt.close()			
		
		print('Winter deviates: ', np.average(winterDists), ' +/- ', np.std(winterDists))
		print('Summer deviates: ', np.average(summerDists), ' +/- ', np.std(summerDists))
		
		
if __name__ == '__main__':
	
	import configHandler as ch
	import sys
	pool = Pool()
	calling_program, masterConigFilename = sys.argv
	
	# Make a config handler object to hold master config file (due to integration with my condor scripts).
	handler = ch.ecumenical(calling_program, masterConigFilename)

	for masterConfig in handler.to_list():
		
		test = mscwFileHolder(masterConfig, pool, masterConfig.general['verbosity'])
		test.extractEventRates()
		test.iterateFit(plot=masterConfig.general['savePlots'], show=masterConfig.general['showPlots'], markAll=False)
		test.excludeRuns()
		test.writeTimeMask()
		test.binRuns()
		test.intraRunMask2(plot=masterConfig.general['savePlots'], show=masterConfig.general['showPlots'])
		test.writeTimeMask()
		#test.plotTimeMasks(show=False)
		#test.plotHistos()
		