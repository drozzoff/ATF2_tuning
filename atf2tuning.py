#!/usr/bin/env python2.7

import numpy as np
import scipy.optimize

from chainmap import ChainMap
from functools import wraps
import pandas as pd
import json
import warnings

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

def gaussian_fit(knob_range, observable_values, **kwargs):
	'''
		Gaussian fit (of the modulation, by default)

		optional parameters:
			init_guess			- [double, double, double]; the list of the initial guesses for the Gaussian (default = [max(observable_values), (knob_range[-1] - knob_range[0]) / 2.0, 0.0])

		Output:
			[
				fit_center,		- double; center of the Gaussian fit
				fit_func		- func; resulting fit function 
			]

	'''
	try: 
		fit_params, pcov = scipy.optimize.curve_fit(gauss, knob_range, observable_values, kwargs.get("init_guess", [max(observable_values), (knob_range[-1] - knob_range[0]) / 2.0, 0.0]), None, True, True, [-np.inf, np.inf], 'trf', gauss_jac)
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

def parabola_fit(knob_range, observable_values, **kwargs):
	'''
		Parabolic fit (of the beamsize, by default)

		optional parameters:
			init_guess				- [double, double, double]; the list of the initial guesses for the Parabola (default = [1e-13, 10, 0.0])
			ignore_restrictions		- bool; if True, the initial fit result check is ignored (default = False) 
											

		Output:
			[
				fit_center,			- double; center of the Gaussian fit
				fit_func			- func; resulting fit function 
			]
	'''
	try: 
		fit_params, pcov = scipy.optimize.curve_fit(parabola, knob_range, observable_values, kwargs.get("init_guess", [1e-13, 10, 0.0]), None, True, True, [-np.inf, np.inf], 'lm', parabola_jac)
	except RuntimeError:
		''' fitting failed, setting to zero '''
		print "\tFitting failed!\n\tSetting to zero."
		return [reduce(lambda x, y: x + y, knob_range) / len(knob_range)]

	if kwargs.get('ignore_restrictions', False):
		return [fit_params[2], lambda x: parabola(x, *fit_params)]
	""" Simple check for the beam size fitting"""
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
			Input:
				beamline				- abstract_machine; input machine that is to be tuned
			
			optional parameters:
				fixed_mode				- bool; if True, IPBSM mode is kept fixe in the tuning (default = False) | is ignored if 'mode' is not set
				initial_errors			- string; name of the file with the initial beamline errors
				ipbsm_errors			- bool;	if True IPBSM dynamic errors are taken into account (default = False)
				mode					- string; mode of the IPBSM | is ignored if fixed_mode = False

			Important (extended definition, covers some of the explicit parameters of the AbstractMachine.iterate_knob()):
			calculation_settings = {
				'lattice'			- string; name of the lattice									| inherits from 'beamline'
				'knobs_info'		- string; knobs definition in the plain format, MAD-X readable. | inherits from 'beamline'
				'post_adj'			- list(func); additional action performed after the observable evaluated.													
													This includes measurement errors, convertion to modulation, custom..
				'duplicate'			- bool; if True, the observable is evaluated twice for each knob value
				'initial_errors'	- string; name of the file with the initial beamline errors 	| is equal to the input 'initial_errors' if it is given as an input
				'ipbsm_errors'		- bool; if True IPBSM dynamic errors are taken into account		| is equal to the input 'ipbsm_errors'
			}
		'''
		assert isinstance(beamline, AbstractMachine), "Input should be of \"abstract_machine\" type"

		self.abs_m, self.knobs_preset, self.duplicate_measurement, self.fixed_mode = beamline, {}, False, kwargs.get('fixed_mode', False)
		self.orbit, self.calculation_settings = Orbit(self.abs_m), self.abs_m.calculation_settings
		self.iteration_log = pd.DataFrame(columns = ['knob', 'keyword', 'fitted_value', 'best_obs', 'mode', 'scan_log', 'setup'])
		self.orbit.read_setup()

		#Setting up the calculation settings
		self.calculation_settings.update({
			'knobs_info'	: self.abs_m.save_knobs_for_madx(),
			'ipbsm_errors'	: kwargs.get('ipbsm_errors', False),
			'duplicate'		: kwargs.get('ipbsm_errors', False)
			})

		if 'initial_errors' in kwargs:
			self.calculation_settings['initial_errors'] = kwargs['initial_errors']

		self.mode = kwargs.get('mode', None)

	def __str__(self):
		return str(self.__dict__)

	__repr__ = __str__

	@property
	def mode(self):
		return self._mode
	
	@mode.setter
	def mode(self, value):
		if value == None: 
			self._mode = None
			return
		if value in ["wire", "6.4", "30", "174"]:
			self._mode = value
			self.calculation_settings['fit'] = parabola_fit if self.mode == 'wire' else gaussian_fit
			self.calculation_settings['post_adj'] = [add_measurement_error(self._mode)] if self.calculation_settings['ipbsm_errors'] else []
			if value in ["6.4", "30", "174"]:
				self.calculation_settings['post_adj'] += [convert_to_modulation(self._mode)]

		else:
			raise ValueError(str(value) + ' - incorrect mode')

	def update_mode(self, data = None):
		'''Updates the mode automatically within the iteration routine'''

		if self.mode == None:
			if data != None and data < (1.4e-6)**2:
				self.mode = "6.4"
				print "Switching to 6.4 degree mode.."
				return
			else:
				self.mode = "wire"
				return

		if self.mode == "wire" and data < (1.4e-6)**2:
			self.mode = "6.4"
			print "Switching to 6.4 degree mode.."
			return

		if self.mode == "6.4" and data > 0.75:
			self.mode = "30"
			print "Switching to 30 degree mode.."
			return

		if self.mode == "30" and data > 0.7:
			self.mode = "174"
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
			if (not self.fixed_mode) and res['keyword'] == "vert. sigma":
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
		'''
			Base function for the beam size tuning
		
			Input:
				var					- string; name of the knob to iterate
				obs					- int, list(list(int)); 
										if obs = 0(default) the vertical beam size is evaluated, 
										if obs = 1 the horizontal beam size is evaluated,

			optional parameters:
				additive			- bool; if True, the knob fitted value is added to the current knob value, otherwise is set to it (default = True)
				points				- int; number of points in the scan (default = 9 (for all the modes except '174'), 21 (for '174' mode))

			*not explicitely used here, but is passed further to machine.construct_mad_request() and affect the result:
				sext_off			- bool; if True, sextupoles are switched off (default = True)
				oct_off				- bool; if True, octupoles are switched off (default = True)
			
			*not explicitely used here, but is passed further to AbstractMachine.iterate_knob():
				plot				- func; if defined the scan results are plotted

		'''
		assert obs in range(2), 'Incorrect parameter' 

		res = self.abs_m.iterate_knob(var, obs, **ChainMap(kwargs, {'points': 21 if self.mode == '174' else 9 }, self.calculation_settings))

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

	@_update_calculation_settings
	@_save_iteration
	def oct2_alignment(self, **kwargs):
		'''
			Iterates the knob AY to find the magnetic center of OCT2

			optional parameters:
				offset_range		- [double, double]; range of the OCT2 misalignments (default = [-1e-3, 1e-3])
				offset_points		- int; number of the points in the scan (default = 7)
				variable			- string; name of the varible in MAD-X that is iterated (default = "oct2dx")
				mode				- string; IPBSM mode: 'wire', '6.4', 30', '174' -> if not given, the current state self.mode is used ->
													if self.mode == None, '174' mode is set
				bba_plot			- func(x, y, fit_func); if defined, the BBA result is plotted

			*not explicitely used here, but is passed further to AbstractMachine.iterate_knob():
				plot				- func; if defined the AY knob scan results are plotted
				knob_range			- list - [double, double]; iteration range for the given knob (default = [-2.0, 2.0])
				points				- int; the number of points in an iteration (default = 9)			-
		'''

		self.mode = kwargs.get("mode", '174' if self.mode == None else self.mode)
		
		if self.mode == '6.4': warnings.warn('Octupole alignment needs 30 or 174 degree mode')
		assert self.mode != "wire", 'Octupole alignment not possible with the Wire Scanner'

		variable = kwargs.get("variable", "oct2dx")
#		res = self.abs_m.iterate_knob('kay', 0, **ChainMap(kwargs, {'points': N_points}, self.calculation_settings))  
		left, right = kwargs.get("offset_range", [-1e-3, 1e-3])
		data = {
			'oct_offset': [],
			'waist_shift': []
		}
		
		for offset in np.linspace(left, right, kwargs.get("offset_points", 7)):
			res = self.abs_m.iterate_knob('kay', 0, **ChainMap(kwargs, {
				'fit': gaussian_fit, 
				'custom_command': variable + " = " + variable + " + ("  + str(offset) + ");\nexec, apply_knobs;",	#should think how to not use MAD-X commands here
				'sext_off': False,
				'oct_off': False
				 }, self.calculation_settings))
			print "OCT2 offset\t" + str(offset) + "\nwaist shift\t" + str(res['fitted_value'])
			data['oct_offset'].append(offset)
			data['waist_shift'].append(res['fitted_value'])
		
		fit_res = parabola_fit(data['oct_offset'], data['waist_shift'], init_guess = [0.0, -1e6, 0.0], ignore_restrictions = True)
		print "\tMagnetic center is " + str(fit_res[0])
		if "bba_plot" in kwargs:
			kwargs.get("bba_plot")(data['oct_offset'], data['waist_shift'], fit_res[1])
		self.knobs_preset['oct2dx'] = fit_res[0]
		return {'knob': "oct2dx", 'keyword': 'oct BBA', 'mode': 'BBA'}

if __name__ == "__main__":
	'''testing routine'''
	'''To do updated'''
