#!/usr/bin/env python2.7

import numpy as np
import scipy.optimize

from chainmap import ChainMap
from functools import wraps
import pandas as pd
import json

from ipbsm import IPBSM
from orbit import Orbit
from machine import AbstractMachine, generate_errors_for_madx

def parabola(x, a, b, c):
	return a + b * (x - c)**2

def parabola_jac(x, a, b, c):
	J = np.zeros((len(x), 3))
	for i in range(len(x)):
		J[i][0] = 1
		J[i][1] = (x[i] - c)**2
		J[i][2] = -2 * b * (x[i] - c)

	return J

def gauss(x, A, sigma, mean):
	return A * np.exp(-(x - mean)**2 / (2 * sigma**2))

def gauss_jac(x, *p):
	A, sigma, mean = p
	J = np.zeros((len(x), 3))
	for i in range(len(x)):
		J[i][0] = gauss(x[i], *p) / A
		J[i][1] = (x[i] - mean)**2 / sigma**3 * gauss(x[i], *p)
		J[i][2] = (x[i] - mean) / sigma**2 * gauss(x[i], *p)

	return J

def convert_to_modulation(mode = None):
	''' custom ATF2 function to use modulation at the IPBSM instead of the beam size '''
	_ipbsm = IPBSM()
	return lambda x: _ipbsm.transform_to_modulation(x, mode)

def add_measurement_error(mode = None):
	'''mode is type of theinstrument used for the measurement'''
	_ipbsm = IPBSM()
	return lambda x: _ipbsm.add_error(x, mode)

def modulation_fit(knob_range, observable_values):
	'''Gaussian fit of the modulation'''
	try: 
		fit_params, pcov = scipy.optimize.curve_fit(gauss, knob_range, observable_values, [max(observable_values), (knob_range[-1] - knob_range[0]) / 2.0, 0.0], None, True, True, [-np.inf, np.inf], 'trf', gauss_jac)
	except RuntimeError: 
		'''The fitting failed'''
		print "\tFitting failed!\n\tSetting to zero."
		return [reduce(lambda x, y: x + y, knob_range) / len(knob_range)]
	'''Checking the boundaries'''
	optimal = fit_params[2]
	if optimal > knob_range[-1]: optimal = knob_range[-1]
	if optimal < knob_range[0]: optimal = knob_range[0]

	print "\tFit is correct!\n\tLargest modulation for " + str(optimal)
	return [fit_params[2], lambda x: gauss(x, *fit_params)]	

def beam_size_fit(knob_range, observable_values):
	''' Parabolic fit of the beam size'''
	try: 
		fit_params, pcov = scipy.optimize.curve_fit(parabola, knob_range, observable_values, [1e-13, 10, 0.0], None, True, True, [-np.inf, np.inf], 'lm', parabola_jac)
	except RuntimeError:
		''' fitting failed, setting to zero '''
		print "\tFitting failed!\n\tSetting to zero."
		return [reduce(lambda x, y: x + y, knob_range) / len(knob_range)]
	if fit_params[1] == 0:
		'''	incorrect fit '''
		print "\tFit is not correct!\n\tSetting to zero."
		return [reduce(lambda x, y: x + y, knob_range) / len(knob_range)]

	if fit_params[1] < 0.0:
		print "\tFit is not correct!\n\tSetting to the boundary."
		if fit_params[2] > 0.0:
			return [knob_range[0]]
		else:
			return [knob_range[1]]
	else:
		
		optimal = fit_params[2]
		'''Checking the boundaries'''
		if optimal > knob_range[-1]: optimal = knob_range[-1]
		if optimal < knob_range[0]: optimal = knob_range[0]

		print "\tFit is correct!\n\tSmallest beam size for " + str(optimal)
		return [optimal, lambda x: parabola(x, *fit_params)]

def atf2_generate_errors():

	#All the magnets
	variables_align = ["QD0DX", "QF1DX", "QD2BDX", "QD2ADX", "QF3DX", "QD4ADX", "QD4BDX", "QF5ADX", "QF5BDX", "QD6DX", "QF7DX", "QD8DX", "QF9ADX", "QF9BDX", "QD10ADX", "QD10BDX", "QM11DX", "QM12DX",
	 "QM13DX", "QM14DX", "QM15DX", "QM16DX", "SF6TTDX", "SF5TTDX", "SD4TTDX", "SF1TTDX", "SD0TTDX", "SK1TTDX", "SK2TTDX", "SK3TTDX", "SK4TTDX", "QD0DY", "QF1DY", "QD2BDY", "QD2ADY", "QF3DY",
	 "QD4ADY", "QD4BDY", "QF5ADY", "QF5BDY", "QD6DY", "QF7DY", "QD8DY", "QF9ADY", "QF9BDY", "QD10ADY", "QD10BDY", "QM11DY", "QM12DY", "QM13DY", "QM14DY", "QM15DY", "QM16DY", "SF6TTDY", "SF5TTDY",
	 "SD4TTDY", "SF1TTDY", "SD0TTDY", "SK1TTDY", "SK2TTDY", "SK3TTDY", "SK4TTDY", "oct1dy", "oct1dx", "oct2dy", "oct2dx"]

	variables_tilt = ["QD0TT", "QF1TT", "QD2ATT", "QD2BTT", "QF3TT", "QD4ATT", "QD4BTT", "QF5ATT", "QF5BTT", "QD6TT", "QF7TT", "QD8TT", "QF9ATT", "QF9BTT", "QD10ATT", "QD10BTT", "QM11TT", "QM12TT",
	 "QM13TT", "QM14TT", "QM15TT", "QM16TT", "SF6TT", "SF5TT", "SD4TT", "SF1TT", "SD0TT", "SK1TT", "SK2TT", "SK3TT", "SK4TT", "oct1tt", "oct2tt"]

	#Base alignment errors for ATF2
	align_error = 1e-4	#in m
	tilt_error = 2e-4	#in rad
	strength_error = 1e-3	#relative error
	bba_error = 1e-4

	variables_strength = ["kqd16x", "kqf17x", "kqd18x", "kqf19x", "kqd20x", "kqf21x", "kqm16ff", "kqm15ff", "kqm14ff", "kqm13ff", "kqm12ff", "kqm11ff", "kqd10bff", "kqd10aff", "kqf9bff", "kqf9aff",
	"kqd8ff", "kqd7ff", "kqd6ff", "kqf5bff", "kqf5aff", "kqd4bff", "kqd4aff", "kqf3ff", "kqd2aff", "kqd2bff", "kqf1ff", "kqd0ff", "ksf6ff", "ksf5ff", "ksd4ff", "ksf1ff", "ksd0ff", "ksk1ff",
	"ksk2ff", "ksk3ff", "ksk4ff", "koct1", "koct2"]

	#Design values
	'''
	init_strength = [-5.159160663, 5.159157136, -3.458719835, 3.027517606, -2.557473406, 4.987927844, 0.0, 5.666474426, -8.99925631, 5.302085025, -1.370824991, 0.8308651056, -1.461012, -1.461218,
	1.907722, 1.9077, -3.044646, 2.77178, -3.034778, 1.894596,  1.894562, -1.495508, -1.495248, 2.784657, -1.459783, -1.001061, 1.661885988, -2.857321177, 220.1178738, 25.80785516, 183.3728196,
	-16.59135058, 30.46821609, -0.009356887575, -1.135747168, -0.3246636104, -0.6069404599, -900.9703714, 17685.41384]
	'''

	init_strength = [-5.159160663, 5.159157136, -3.458719835, 3.027517606, -2.557473406, 4.987927844, 0.0, 5.666474426, -8.99925631, 5.302085025, -1.370824991, 0.8308651056, -1.461012, -1.461218,
	1.907722, 1.9077, -3.044646, 2.77178, -3.034778, 1.894596,  1.894562, -1.495508, -1.495248, 2.784657, -1.459783, -1.001061, 1.661885988, -2.857321177, 220.1178738, 25.80785516, 183.3728196,
	-16.59135058, 30.46821609, -0.009356887575, -1.135747168, -0.3246636104, -0.6069404599, 0.0, 7309.9998]

	norm_sext_align_variables = [["SD0TTDX", "SD0TTDY"], ["SF1TTDX", "SF1TTDY"], ["SD4TTDX", "SD4TTDY"], ["SF5TTDX", "SF5TTDY"], ["SF6TTDX", "SF6TTDY"]]
	oct_align_variables = [["oct1dx", "oct1dy"], ["oct1dx", "oct1dy"]]

	variables = [] * (len(variables_align) + len(variables_tilt) + len(variables_strength))
	for x in variables_align:
		variables.append([x, align_error, 3 * align_error, False])
	for x in variables_tilt:
		variables.append([x, tilt_error, tilt_error, False])

	for i in range(len(variables_strength)):
		error = abs(init_strength[i] * strength_error)
		variables.append([variables_strength[i], error, 3 * error, True])
	#	Initializing the setup
	#Generating initial errors

	generate_errors_for_madx(variables)

class ATF2Tuning(object):
	"""
		Class that is used for the tuning. It has a general tuning routine, custom ones, such as wire scan and QD0 roll scan as well as orbit/dispersion measurement/correction
	"""
	def __init__(self, beamline, **kwargs):
		'''
			Constructor input:
				beamline	- abstract_machine; input machine that is to be tuned
		'''
		assert isinstance(beamline, AbstractMachine), "Input should be of \"abstract_machine\" type"

		self.abs_m, self.mode, self.knobs_preset, self.duplicate_measurement = beamline, None, {}, False
		self.orbit, self.calculation_settings = Orbit(self.abs_m), self.abs_m.calculation_settings, 
		self.iteration_log = pd.DataFrame(columns = ['knob', 'keyword', 'fitted_value', 'best_obs', 'mode', 'scan_log', 'setup'])
		self.orbit.read_setup()

		#Setting up the calculation settings
		self.calculation_settings['knobs_info'] = self.abs_m.save_knobs_for_madx()
		if kwargs.get('ipbsm_errors', False): 
			self.calculation_settings['post_adj'] = [add_measurement_error(self.mode)]
			self.calculation_settings['duplicate'] = True
		else:
			self.calculation_settings['post_adj'] = []

		if 'initial_errors' in kwargs:
			self.calculation_settings['initial_errors'] = kwargs['initial_errors']

		self.update_mode()

	def __str__(self):
		return str(self.__dict__)

	__repr__ = __str__

	def update_mode(self, data = None):
		'''Updates the mode automatically within the iteration routine'''
#		modes_list = ["wire", "8", "30", "174"]
		modes_list = ["wire", "6.4", "30", "174"]
		if self.mode == None:
			if data == None:
				self.mode = "wire"
				self.calculation_settings['fit'] = beam_size_fit
				return
			elif data < (1.4e-6)**2:
#				self.mode = "8"
				self.mode = "6.4"	#Using 6.4 degree mode as a default initial stage of the IP-BSM tuning
				self.calculation_settings['post_adj'] += [convert_to_modulation(self.mode)]
				self.calculation_settings['fit'] = modulation_fit
				return
			else:
				self.mode = "wire"
				return

		if self.mode == "wire" and data < (1.4e-6)**2:
#			self.mode = "8"
#			print "Switching to 8 degree mode.."
			self.mode = "6.4"
			self.calculation_settings['post_adj'] += [convert_to_modulation(self.mode)]
			self.calculation_settings['fit'] = modulation_fit
			print "Switching to 6.4 degree mode.."
			return
		if self.mode == "6.4" and data > 0.75:
			self.mode = "30"
			self.calculation_settings['post_adj'][-1] = convert_to_modulation(self.mode)
			print "Switching to 30 degree mode.."
			return
		if self.mode == "30" and data > 0.7:
			self.mode = "174"
			self.calculation_settings['post_adj'][-1] = convert_to_modulation(self.mode)
			print "Switching to 174 degree mode.."
			return
		return

	def _update_calculation_settings(func):
		@wraps(func)
		def wrapper(self, *arg, **kwargs):
			res = func(self, *arg, **kwargs)
			self.calculation_settings['knobs_preset'] = self.knobs_preset
			return res
		return wrapper

	def _update_mode(func):
		@wraps(func)
		def wrapper(self, *arg, **kwargs):
			res = func(self, *arg, **kwargs)
			if res['keyword'] == "vert. sigma":
				self.update_mode(res['best_obs'])
			return res
		return wrapper

	def _save_iteration(func):
		@wraps(func)
		def wrapper(self, *arg, **kwargs):
			res, setup = func(self, *arg, **kwargs), {}
			try:
				setup = self.calculation_settings['knobs_preset']
			except KeyError: pass

			self.iteration_log = self.iteration_log.append(dict(ChainMap(res, {'setup': json.dumps({'setup': setup})} )), ignore_index = True)
			print self.iteration_log
			return res
		return wrapper

	@_update_calculation_settings
	@_update_mode
	@_save_iteration
	def tune(self, var, obs = 0, **kwargs):
		'''Base function for the beam size tuning'''
		assert obs in range(2), 'Incorrect parameter' 
		N_points = 9
		if self.mode == "174": N_points = 21

		res = self.abs_m.iterate_knob(var, obs, **ChainMap(kwargs, {'points': N_points}, self.calculation_settings)) # [knob_amp, best_value]

		if (var in self.knobs_preset) and kwargs.get('additive', True):
			self.knobs_preset[var] += res['fitted_value']
		else:
			self.knobs_preset[var] = res['fitted_value']

		return dict(ChainMap({'knob': var, 'keyword': 'hor. sigma' if obs else 'vert. sigma', 'mode': self.mode}, res))

	def measure_orbit(self, **kwargs):
		'''Orbit measurement | Input is the calculation settings dict'''
		self.orbit.measure_orbit(**ChainMap(kwargs, self.calculation_settings))

	@_update_calculation_settings
	@_save_iteration
	def correct_orbit(self, **kwargs):
		'''Orbit correction | keeping the correctors setup in memory'''
		self.knobs_preset = dict(ChainMap(self.orbit.correct_orbit(**ChainMap(kwargs, self.calculation_settings)), self.knobs_preset))
		return {'knob': "orbit corr", 'keyword': 'hor. and vert. BPMs', 'mode': 'orbit'}

	"""Files handling"""
	def save_log(self, filename = "log.pkl"):
		self.iteration_log.to_pickle(filename)

	def save_setup(self, filename = "setup.json"):
		with open(filename, 'w') as f:
			json.dump(self.knobs_preset, f)

	def read_setup(self, filename = "setup.json"):
		with open(filename, 'r') as f:
			self.knobs_preset = json.load(f)
		self.calculation_settings['knobs_preset'] = self.knobs_preset


	'''07.02.2021

		Octupole alignment is not yet updated
	'''
	def oct2_align(self, variable, **kwargs):
		'''alignment routine of OCT2 based on the waist shift at the IP'''
		default_value = {'step': 21, 'shift_range': [-1e-3, 1e-3], 'plot_option': False, 'add_variable': None, 'add_shift': None}
		adjust = lambda x: kwargs.get(x, default_value[x])
		for x in default_value:
			default_value[x] = adjust(x)

		best_value, data_log = [], []
		N_points, default_strategy, plot_option = 7, 3, True
		
#		if self.mode == "174": N_points = 15

		def update_location(value):
			with open("bba_alignment.madx", 'w') as f:
				print >> f, variable + " = " + variable + " + ("  + str(value) + ");"
				if default_value['add_variable'] == None: return
				print >> f, default_value['add_variable']  + " = " + default_value['add_variable'] + " + ("  + str(default_value['add_shift']) + ");"
		peaks = []
		shift_range = default_value['shift_range']
#		print shift_range, default_value['step']
		for shift in list(np.linspace(shift_range[0], shift_range[1], default_value['step'])):
#		for shift in [-0.000125]:
			update_location(shift)
			knob_optimum, __, __ = self.abs_m.iterate_knob('kay_2', self.observable, default_value['plot_option'], False, [-0.2, 0.2], default_strategy, self.mode, N_points, False)
			peaks.append(knob_optimum)
			print shift, knob_optimum

		print peaks

if __name__ == "__main__":
	'''testing routine'''
	'''To do updated'''
