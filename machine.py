#!/usr/bin/env python2.7

import numpy as np

from coordinatecore import CoordinateCore

from numpy import linalg as LA

import random
import json
import subprocess
from functools import wraps


#Python 3.2+
#from collections import ChainMap

#Python 2.7 | custom chainmap module
from chainmap import ChainMap

def write_mad_command(*args, **kwargs):
	filename = kwargs.get('filename', "temp.madx")
	def inner_func(func):
		@wraps(func)
		def wrapper(*arg, **kwargs):
			command = func(*arg, **kwargs)
			with open(filename, 'w') as f:
				print >> f, command
			return command
		return wrapper
	return inner_func


def mad_execute(flag = True):
	'''Gotta be modified to use PyMAD'''
	if flag:
		subprocess.call("./madx < job.madx > trash.out", shell = True)
	else:
		subprocess.call("./madx < job.madx", shell = True)


@write_mad_command(filename = "initial_errors.madx")
def generate_errors_for_madx(variables_list, append = False):
	'''
	Generating the errors for MadX
	variables_list is a list with pairs [variable, sigma, boundary, modify init value] - name of the variable and the rms errors

	'''
	s = ""
	for x in variables_list:
		dd = 1.1 * x[2]
		while abs(dd) > x[2]:	#rule taken from the older codes
			dd = random.gauss(0, x[1])
		if x[3]:
			s += x[0] + " = " + x[0] + " + (" + str(dd) + ");\n"
		else:
			s += x[0] + " = " + str(dd) + ";\n"
	return s


class Knob(object):
	"""
		Structure that is used to store the knob data
	"""
	def __init__(self, **kwargs):
		'''
		Constructor input:
			name		- string; name of the vector
			vector		- list(float); vector of the associated deviations (shifts, strengths deviation..)
			variables	- list(string); vector of the associated variables (names of the corresponding shifts, strength deviation..)
		'''
		
		if 'filename' in kwargs:
			self.read(kwargs.get("filename"))
			return

		self.name, self.variables, self.vector, self.label = kwargs.get('name', "noname"), kwargs.get('variables', []), kwargs.get('vector', []), kwargs.get('label', '')

	def __str__(self):
		return json.dumps(self.__dict__)

	def __mul__(self, factor):
		'''Scales the knob.vector'''
		self.vector = [x * factor for x in self.vector]

	def save(self):
		'''Saves the knob parameters to a file with a name based on the name of the knob'''
		with open(str(self.name) + ".knob", 'w') as f:
			json.dump(self.__dict__, f)
		return None

	def read(self, filename):
		with open(filename, "r") as f:
			self.__dict__ = json.load(f)

	__repr__ = __str__

class Machine(CoordinateCore):
	"""
	Main class for the knobs construction/knobs iteration
	"""
	def __init__(self, sigmaInitial, gaussdpp, **kwargs):
		'''
			Constructor input
				sigmaInitial	- list; list of the sigmax at the beam line entrance
				gaussdpp		- bool; If True, Gaussian distribution for delta is applied, otherwise - Flat-top
			
			Optional input:
				order			- int; The order of the DA map used in the calculations (default = 2)
				method			- int; If 1 - Mapclass is used to evaluate the distribution, if 0 - internal routines are used (are slower) (default = 1)
				name			- string; The name of the beamline (default= default)
			
		'''

#		self.map_file = map_file
		self.order, self.method, self.knobs, self.name= kwargs.get('order', 2), kwargs.get('method', 0), None, kwargs.get('name', "default")
		self.construct_mad_request(lattice = self.name)
		super(Machine, self).__init__(sigmaInitial, gaussdpp, **kwargs)
#		self.tuning_order, self.Nparticles, self.method = kwargs.get('tuning_order', 2), , kwargs.get('method', 1)
		
#		exit()
		self.construct_mad_request(lattice = self.name)

		#variables needed for the knobs construction/tuning
		self.variables, self.weights = None, None

	def __str__(self):
		return str(self.__dict__)

	__repr__ = __str__

	def execute_mad_and_update_map(func):
		@wraps(func)
		def wrapper(self, *args, **kwargs):
			func(self, *args, **kwargs)
			mad_execute()
			self.update()
		return wrapper

	@execute_mad_and_update_map
	@write_mad_command(filename = 'input.madx')
	def construct_mad_request(self, **kwargs):
		'''
			Writing the commands to input.madx file based on the priority
			1. initial errors
			2. knob file
			3. initial knob values
			4. rest of the initial changes
			5. apply the knob
		
		'''
##########################################################################################
		'''initialization of the sequence'''
		command = 'use, period = ' + kwargs.get('lattice', 'default') + ';\n'
##########################################################################################
		'''Setting up the machine'''

		''''Initial errors of the lattice'''
		if 'initial_errors' in kwargs:
			command += 'call, file = ' + kwargs['initial_errors'] + ';\n'
			'''bba included'''
			if 'bba_settings' in kwargs:
				command += kwargs['bba_settings']
			command += 'exec, apply_errors;\n'
		
		'''Nonlinear elements off for the tuning'''
		if kwargs.get('sext_off', True):
			command += 'exec, switch_sextupoles_off;\n'
		if kwargs.get('oct_off', True):
			command += 'exec, switch_octupoles_off;\n'

		'''The knobs definition and current values | +corrector settings'''
		if 'knobs_info' in kwargs:
			command += kwargs['knobs_info']
		if 'knobs_preset' in kwargs:
			knob_list = kwargs['knobs_preset']
			for app_knob in knob_list:
				command += app_knob + ' = ' + str(knob_list[app_knob]) + ';\n'
			command += 'exec, apply_knobs;\n'
			command += 'exec, reset_knobs;\n'
##########################################################################################
		'''Separately applying a knob | could be merged with the upper sections'''
		if 'apply_knob' in kwargs:	#apply_knob is expected to be a dictionary
			knob_list = kwargs['apply_knob']
			for app_knob in knob_list:
				command += app_knob + ' = ' + str(knob_list[app_knob]) + ';\n'

			command += "exec, apply_knobs;\n"
##########################################################################################
		'''Custom command to be executed in MadX | Mainly used when constructing the knobs'''
		if 'custom_command' in kwargs:
				command += kwargs['custom_command'] + '\n' 
##########################################################################################
		'''orbit handling using MadX built in tools'''
		if kwargs.get('correct_orbit', False):
			command += "exec, twiss_macro;\n"
			command += "exec, orbit_correction;\n"

		if kwargs.get('measure_orbit', False):
			command += 'exec, twiss_macro;\n'
			command += 'exec, write_bpm_data;\n'
##########################################################################################
		return command

	def construct_sigma_matrix(self):
		'''
			Returns the calculated Sigma matrix
		'''
		Sigma_m = np.zeros((5, 5))
		for i in range(5):
			for j in range(5):
				Sigma_m[i][j] = self.get_sigma([i, j], self.method)

		return Sigma_m

	def precalculate_terms(self, list_of_terms):
		'''
			Calculates (machine.init_terms) and returns the sigma for a given coordinates combinations. It is used when the sigma is calculated for more than 3 coordinates (otherwise the sigma matrix is used).
			It is mainly used for nonlinear knobs construction
			
			Input:
				list_of_terms	- list(list(int)); holds the list of terms vectors
					Ex: list_of_terms = [[0,1], [2,2]] means that we want sigma_x_x' and sigma_y_y
		'''

		#Updating the map in the memory
		self.construct_mad_request(lattice = self.name)

		self.init_terms = list(map(lambda x: self.get_sigma(x, self.method), list_of_terms))
		return self.init_terms


	def build_response_matrix(self, deviation, **kwargs):
		'''
			Calculates and returns the response of the terms in list_of_terms (or in Sigma-matrix) on change of the values of parameters in machine.variables (has to be defined before the execution) 
			
			The default response matrix is build in the form:
			[[d/dx_1(sigma_{x,x}),	d/dx_2(sigma_{x,x}),	d/dx_3(sigma_{x,x})...]
			 [d/dx_1(sigma_{x,x'}),	d/dx_2(sigma_{x,x'}),	d/dx_3(sigma_{x,x'})...]
			 ...
			 [d/dx_1(sigma_{delta,delta}),	d/dx_2(sigma_{delta,delta}),	d/dx_3(sigma_{delta,delta})...]]
			
			Input:
				list_of_terms	- list(list(float); list of the terms to calcuate he response of
					Ex: list_of_terms = [[0,1], [2,2]] means that we want sigma_x_x' and sigma_y_y to be recorded to the response matrix
				weight_type		- int; there are several predefined weights for the response matrix calculation
					1.(default)		R_ij = d/dx sigma_i_j / sqrt(sigma_i_i * sigma_j_j) 
					2.				R_ij = d/dx sigma_i_j / sqrt(sigma_j_j) 
					3.				R_ij = d/dx sigma_i_j / sigma_j_j 

					*The function ignores weight_type flag if the weights array is defined

			!!In the case when 3 variables sigma is needed, the precalculate_terms() should be called before using this function to calculate the initial values.

			17.09.2021: Description in outdated
		'''

		#system("cp knobs_contruction.madx input.madx")
		
		list_of_terms, weight_type, relative = kwargs.get('list_of_terms', None), kwargs.get('weight_type', 1),  kwargs.get('relative', True),
		self.variables = kwargs.get('variables', self.variables)

		assert self.variables != None, "Knobs construction parameters are not set!"

		self.construct_mad_request(lattice = self.name, **kwargs)

		self.sigma_m_0 = self.construct_sigma_matrix()

		print self.sigma_m_0

		if list_of_terms == None:
			self.R = np.zeros((25,len(self.variables)))
		else:
			self.R = np.zeros((len(list_of_terms), len(self.variables)))

		x = 0
		for var in self.variables:
			print "Changing ", var, " -  Sigma matrix is:"
			#delta = self.deviation
			command = None
			if relative == True:
				command = var + "=" + var + " * (1 +" + str(deviation) + ");\n"
			else:
				command =  var + "=" + str(deviation) + ";\n"
			command += 'exec, apply_knobs;\n'

			self.construct_mad_request(lattice = self.name, custom_command = command, **kwargs)

			sigma_m = self.construct_sigma_matrix()
			print sigma_m

			if list_of_terms == None:
				#default R-matrix construction 
				for i in range(5):
					for j in range(5):
						#this wont work for the relative shift
						if self.weights == None:
							self.R[5 * i + j][x] = (sigma_m[j][i] - self.sigma_m_0[j][i]) / deviation * (1 / self.sigma_m_0[j][i])
						else:
							self.R[5 * i + j][x] = (sigma_m[j][i] - self.sigma_m_0[j][i]) / deviation * self.weights[5 * i + j]
				x += 1
			else:
				#Constructing the R-matrix from the defined variables
				for i in range(len(list_of_terms)):
					#default weight in 2d case is 1 / sqrt(sigma_j_j * sigma_i_i)
					a = list_of_terms[i][0]
					b = list_of_terms[i][1]
					if len(list_of_terms[0]) == 2:
						if self.weights == None:
							if weight_type == 1:
								self.R[i][x] = (sigma_m[a][b] - self.sigma_m_0[a][b]) / np.sqrt(sigma_m[a][a] * sigma_m[b][b]) / deviation
							elif weight_type == 2:
								self.R[i][x] = (sigma_m[a][b] - self.sigma_m_0[a][b]) / np.sqrt(sigma_m[b][b]) / deviation
							elif weight_type == 3:
								self.R[i][x] = (sigma_m[a][b] - self.sigma_m_0[a][b]) / sigma_m[b][b] / deviation
						else:
							self.R[i][x] = (sigma_m[a][b] - self.sigma_m_0[a][b]) * self.weights[i] / deviation
					elif len(list_of_terms[i]) == 3:
						c = list_of_terms[i][2]
						term0 = self.init_terms[i]
						term = self.get_sigma(list_of_terms[i], self.mathod)
						if self.weights == None:
							if weight_type == 1:
								self.R[i][x] = (term - term0) / np.sqrt(sigma_m[a][a] * sigma_m[b][b] * sigma_m[c][c]) / deviation
							elif weight_type == 2:
								self.R[i][x] = (term - term0) / np.sqrt(sigma_m[b][b] * sigma_m[c][c]) / deviation
							elif weight_type == 3:
								self.R[i][x] = (term - term0) / sigma_m[b][b] / sigma_m[c][c] / deviation
						else:
							self.R[i][x] = (term - term_0) * self.weights[i] / deviation
					else:
						raise ValueError("Unsupported order")
				
				x += 1
		
		print 'Response matrix is:\n', self.R
		
		return self.R

	def construct_knobs(self, knobs_name, label, **kwargs):
		'''
			Constructs and returns the knobs.
			
			Notes:
				It uses already precalculated response matrix <-- should include an Exception
				The length of knobs_name should correpond to the list_of_terms in the response matrix construction routine

			Input:
				knobs_name		- list(string); names to call the constructed knobs
									len(knobs_name) has to be the same as len(self.list_of_terms) <-- should include an Exception
		'''

		assert 'R' in self.__dict__, "Response Matrix is missing, run ``build_response_matrix()'' before"

		U, E, V_t = LA.svd(self.R)

		V = V_t.transpose()
		
		print "Matrix V_t is\n", V_t
		print "Matrix U is\n", U

		if kwargs.get('overwrite_knobs', False) or (self.knobs == None):
			self.knobs = []

		for i in range(len(self.R)):
			self.knobs.append(Knob(variables = self.variables, vector = list(V_t[i]), name = knobs_name[i], label = label))

		return self.knobs

	def save_knobs(self):
		''' Writes the knobs in the memory to a file '''
		__ = map(lambda x: x.save(), self.knobs)

	@write_mad_command(filename = 'knob.madx')
	def save_knobs_for_madx(self):
		'''
			Saves the knobs to a file, which can be called from MadX || Gotta be rewritten in the future to contain a string only

			Input:
				fileName	- string; name of the file
		'''
		#iterating over the labels associated with the knobs
		if self.knobs == None: return ""

		command, labels = "", []

		for x in self.knobs:
			if not (x.label in labels): labels.append(x.label)
		
		for l in labels:
			#the knobs with the same label are assumed to be constructed on the same variables
			current_knobs = list(filter(lambda x: x.label == l, self.knobs))	
			variables = current_knobs[0].variables

			for i in range(len(variables)):
				command += variables[i] + ' := '
				for j in range(len(current_knobs)):
					vector = current_knobs[j].vector
					name = current_knobs[j].name
					if j == 0:
						command +=  name + ' * (' + str(vector[i]) +  ')'
					else:
						command += ' + ' + name + ' * (' + str(vector[i]) + ')'
				command += ';\n'
		return command

class AbstractMachine(Machine):
	"""
		Extension of the machine class with knobs handling routines
	"""
	def __init__(self, sigmaInitial, gaussdpp, **kwargs):
		'''
			!!!have to be updated -> check base class
			Constructor input
				tuning_order	- int; maximum order of the map to be used in the beam size calculation
				Nparticles		- int; number of particles to be used in the tracking if the internal tracking of coordinates_core is used
				method			- int(0 or 1); sigma calculation method (coordinates_core)

		'''
		super(AbstractMachine, self).__init__(sigmaInitial, gaussdpp, **kwargs)
		self.beam_size_errors, self.calculation_settings = kwargs.get('beam_size_errors', False), {}

		self.calculation_settings['lattice'] = self.name

		if self.knobs != None:
			self.calculation_settings['knobs_info'] = self.save_knobs_for_madx()

	def __str__(self):
		return str(self.__dict__)

	__repr__ = __str__
	"""Knob adjustments"""
	def knob_check(self, knob_name, **kwargs):
		'''
			Iterates the knob and checks the change of the observables

			Input:
				knob		- string; name of the knob applied
				observables	- list(list(int)); list of the values to be measured for the given value of the knob 
					(prints the difference from the initial values)
					Ex: observables = [[0,1], [2,2]] means that we want sigma_x_x' and sigma_y_y contributions
				range_scale	- float; factor to modify the default iteration range [-1.0, 1.0]
				fileName	- string; name of the file to save knob iteration information
		'''
		assert knob_name in list(map(lambda x: x.name, self.knobs)), 'knob_check() - Knob "' + knob_name + '" does not exist'

		observables, range_scale, fileName = kwargs.get('observables', None), kwargs.get('range_scale', 1.0), kwargs.get('filename', None)
		knob_files = kwargs.get('knob_fiels', [])

		self.construct_mad_request(lattice = self.name)

		self.sigma_m_0 = self.construct_sigma_matrix()
		
		knob_range = np.array(np.linspace(-1e0, +1e0, kwargs.get('points', 21))) * range_scale
		
		if len(observables[0]) == 3: self.precalculate_terms(observables)
		
		log_data = knob_name
		for x in observables: log_data += '\t' + str(x)
		print log_data
		log_data += '\n'

		for value in knob_range:
			ln = str(value)

			self.construct_mad_request(apply_knob = {knob_name: value}, **ChainMap(kwargs, self.calculation_settings))

			sigma_m = self.construct_sigma_matrix()

			for i in range(len(observables)):
				a, b, x = observables[i][0], observables[i][1], 0.0

				if len(observables[i]) == 2: 
					x = (sigma_m[a][b] - self.sigma_m_0[a][b]) / np.sqrt(sigma_m[b][b])
				elif len(observables[i]) == 3:
					c = observables[i][2]
					x = (self.get_sigma(observables[i], self.method) - self.init_terms[i]) / np.sqrt(sigma_m[b][b] * sigma_m[c][c])
				else:
					raise ValueError("Unsupported dimension")
				ln += '\t' + str(x)
			print ln
			log_data += ln + '\n'

		return log_data


	def normalize_knob(self, knob_name, observable, **kwargs):
		'''
			Modifies the knob in a such a way that 1.0 knob amplitude would give target value for the observable
			The knob is applied once at the amplitude 'value' 

			
			Input:
				knob_name	- string; name of the knob to be normalized
				observable	- list(int); main covariance to be observed
				value		- float; knob amplitude
				target		- float; target beam size contribution in [nm]

			By default, normalizing the knobs in way that amplitude of 1 means 100 nm contribution.
			Using the value of the observable at the amplitude of 1 * value_modifier

		'''

		amplitude, target = kwargs.get('amplitude', 1.0),  kwargs.get('target', 1.0)
		
		assert knob_name in list(map(lambda x: x.name, self.knobs)), 'normalize_knob() - Knob "' + knob_name + '" does not exist'

		self.construct_mad_request(lattice = self.name)
		
		in_term = 0.0
		if len(observable) == 3:
			in_term = self.get_sigma(observable, self.mathod)
		else:
			self.sigma_m_0 = self.construct_sigma_matrix()

		log_data = "Applying the knob: " + knob_name + ':' + str(amplitude)

		self.construct_mad_request(apply_knob = {knob_name: amplitude}, **ChainMap(kwargs, self.calculation_settings))

		a, b, correction, sigma_m = observable[0], observable[1], 0.0, self.construct_sigma_matrix()
	
		if len(observable) == 2: 
			correction = np.abs((sigma_m[a][b] - self.sigma_m_0[a][b]) / np.sqrt(sigma_m[b][b]))
		elif len(observable) == 3:
			c = observable[2]
			correction = np.abs((self.get_sigma(observable, self.method) - in_term) / np.sqrt(sigma_m[b][b] * sigma_m[c][c]))
		else:
			pass
			raise Exception("normalize_knob() - Unsupported dimension (yet?)")
		
		for i in self.knobs:
			if i.name == knob_name:
				log_data +=  " | Correction factor is " + str(target / correction * amplitude)
				i.vector = list(map(lambda x: x * target / correction * amplitude, i.vector))
		print log_data
#		return log_data
		
	def _get_obs(self, observable, **kwargs):
		'''Returns sigma or some specific terms'''
		self.construct_mad_request(**ChainMap(kwargs, self.calculation_settings))	
		return self.get_obs(observable = observable)

	def iterate_knob(self, knob_name, observable, **kwargs):
		'''
			Iterates the knob in the range (default is [-2.0, 2.0]) and checking the vertical beam size. After matches the knob value to get the smallest beam size.
			
			Input parameters:
				knob_name			- string; name of the knob to iterate
				observable			- variable to be checked for the iteration

				show_plot			- bool; plotting options
				sext_flag			- bool; if True, normal sextupoles are ON, False - OFF
				range				- list(float); iteration range
				strategy			- int(1, 2 or 3); different strategies for the observable (beam_size) tuning:
										1: Beam size tuning | Parabola fit
										2: Same as (1) but the beam size is measured twice for each amplitude | Parabola fit
										3: Modulation is used (2 measurements for each amplitude) | Gaussian fit
				mode				- string; definition of the measurement mode, such as "wire", "8", "30" and "174" degree modes
				Number_of_points	- int;  number of points used within the given range

			***
			If beam_size_errors is True, the strategies, where 2 measurements are expected will be ignored
		'''
		
		post_adj, knob_range = kwargs.get('post_adj', []), kwargs.get('knob_range', [-2.0, 2.0])
		knob_range = np.array(np.linspace(knob_range[0], knob_range[1], kwargs.get('points', 9)))

		if kwargs.get('duplicate', False):
			knob_range = np.sort(np.concatenate((knob_range, knob_range)))

		observable_values = []

		print "Iterating the knob " + knob_name + ":"
		if observable == 0:
			print "\tMeasuring the vertical beam size"
		if observable == 1:
			print "\tMeasuring the horizontal beam size"
		
		for x in knob_range:
			self.construct_mad_request(apply_knob = {knob_name: x}, **ChainMap(kwargs, self.calculation_settings))
			
			obs = self.get_obs(observable = observable)
			print "\t" ,  x, "\t", obs,
			for f in post_adj:
				obs = f(obs)
				print "->\t", obs,
			observable_values.append(obs)
			print ""

		fit_result = None	
		if 'fit' in kwargs: fit_result = kwargs['fit'](knob_range, observable_values)
		if ('plot' in kwargs):
			if (fit_result != None) and len(fit_result) == 2:
				kwargs['plot'](knob_range, observable_values, fit_result[1])
			else:
				kwargs['plot'](knob_range, observable_values)

		iter_data = json.dumps({'knob_range': knob_range.tolist(), 'obs_data': observable_values})

		if fit_result != None:
			if len(fit_result) == 2:
				return {'scan_log': iter_data, 'fitted_value': fit_result[0], 'best_obs': fit_result[1](fit_result[0])}
			else:
				'''evaluating the obs'''
				self.construct_mad_request(apply_knob = {knob_name: fit_result[0]}, **ChainMap(kwargs, self.calculation_settings))
				return {'scan_log': iter_data, 'fitted_value': fit_result[0], 'best_obs': self.get_obs(observable = observable)}
		else: return {'scan_log': iter_data}


'''Testing routine'''
if __name__ == "__main__":
	'''Testing routine'''
	'''To be updated'''