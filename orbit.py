#!/usr/bin/env python2.7

import sys
sys.path.append('/afs/cern.ch/eng/sl/lintrack/MAPCLASS/MapClass2/')
from metaclass2 import twiss2

import numpy as np
import copy

import json
from machine import AbstractMachine

#Python 3.2+
#from collections import ChainMap

#Python 2.7 | custom chainmap module
from chainmap import ChainMap

class Orbit(object):
	'''Base class to work with the orbit (BPMs and corrs)'''
	def __init__(self, beamline, **kwargs):
		''' '''
		assert isinstance(beamline, AbstractMachine), "Input should be of \"machine\" type"

		self.beamline, self.calculation_settings = beamline, beamline.calculation_settings
		if 'initial_errors' in kwargs:
			self.calculation_settings['initial_errors'] = kwargs['initial_errors']

	def __str__(self):
		return str(self.__dict__)

	__repr__ = __str__

	def read_bpms_setup(self, filename = "setup.bpms"):
		with open(filename, 'r') as f:
			self.bpm_list = json.load(f)

	def read_kickers_setup(self, filename = "setup.corrs"):
		with open(filename, 'r') as f:
			self.corr_list = json.load(f)

	def write_bpms_setup(self, filename = "setup.bpms"):
		with open(filename, 'w') as f:
			json.dump(self.bpm_list, f)

	def write_kickers_setup(self, filename = "setup.corrs"):
		with open(filename, 'w') as f:
			json.dump(self.corr_list, f)
	
	def read_setup(self, **kwargs):
		''' '''
		self.read_bpms_setup(kwargs.get('bpm_files', "setup.bpms"))
		self.read_kickers_setup(kwargs.get('corr_files', "setup.corrs"))

		'''extrackting bpm locations'''
		for element in twiss2("bds.twiss").elems:
			if element.NAME in self.bpm_list:
				self.bpm_list[element.NAME]['s'] = element.S
			if element.NAME.lower() in self.bpm_list:
				self.bpm_list[element.NAME.lower()]['s'] = element.S

			if element.NAME in self.corr_list:
				self.corr_list[element.NAME]['s'] = element.S
			if element.NAME.lower() in self.corr_list:
				self.corr_list[element.NAME.lower()]['s'] = element.S

	def read_orbit(self, **kwargs):
		'''reading the correctots strengths and orbit from the files'''
		for bpm in twiss2(kwargs.get('bpm_files', ["mon_x.out", "mon_y.out"])[0]).elems:
			if bpm.NAME in self.bpm_list:
				self.bpm_list[bpm.NAME]['x'] = bpm.X
		for bpm in twiss2(kwargs.get('bpm_files', ["mon_x.out", "mon_y.out"])[1]).elems:
			if bpm.NAME in self.bpm_list:
				self.bpm_list[bpm.NAME]['y'] = bpm.Y


		for corr in twiss2(kwargs.get('corr_files', ["corr_x.out", "corr_y.out"])[0]).elems:
			if corr.NAME in self.corr_list:
				self.corr_list[corr.NAME]['hkick'] = corr['PX.CORRECTION']
		for corr in twiss2(kwargs.get('corr_files', ["corr_x.out", "corr_y.out"])[1]).elems:
			if corr.NAME in self.corr_list:
				self.corr_list[corr.NAME]['vkick'] = corr['PY.CORRECTION']

	def measure_orbit(self, **kwargs):
		'''Measuring the current orbit.'''
		self.beamline.construct_mad_request(measure_orbit = True, **ChainMap(kwargs, self.calculation_settings))

		self.read_orbit(**kwargs)

		if 'plot' in kwargs:
			kwargs['plot'](self.bpm_list)

	def correct_orbit(self, **kwargs):
		'''Correcting the orbit and measuring it'''
		self.beamline.construct_mad_request(correct_orbit = True, **ChainMap(kwargs, self.calculation_settings))
		self.read_orbit(**kwargs)

		if 'plot' in kwargs:
			kwargs['plot'](self.bpm_list)
		return self.to_knobs_preset()

	def to_knobs_preset(self):
		''' '''
		result = {}
		for corr in self.corr_list:
			if self.corr_list[corr]['status'] == "on":
				kick = 0.0
				if (self.corr_list[corr]['plane'] == "hor") and ('hkick' in self.corr_list[corr]):
					kick = self.corr_list[corr]['hkick']
				if (self.corr_list[corr]['plane'] == "vert") and ('vkick' in self.corr_list[corr]):
					kick = self.corr_list[corr]['vkick']
				
				result[corr + "->kick"] = kick
		return result

if __name__ == "__main__":
	'''Testing routine'''
	'''To be updated'''