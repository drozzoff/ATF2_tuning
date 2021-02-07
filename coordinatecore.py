#!/usr/bin/env python2.7

import sys
#path to Mapclass libs on afs
sys.path.append('/afs/cern.ch/eng/sl/lintrack/MAPCLASS/MapClass2/')	

import numpy as np

import mapclass

class CoordinateCore(object):
	"""

		Is used to calculate the coordinates of an individual particle or the total sigmas after passing through the lattice.
		2 approaches are implemented:
			With MAPCLASS - works well for sigma of 3 or less variables.
			With MAPCLASS + custom tracker - should be used for sigma of 4 or more variables (Octupole knobs...). Map of the maximum order of 3 is accepted in this case.
				This method has a drawback -> the orbit deviation is ignored. When it is an important, the beam size obtained is not correct.
				**Should be corrected**
			

		Internally uses the variable sigma_initial, which has to be defined externally in the program
		
		Uses "fort.18" file to read the map
	"""
	state = ['x', 'px', 'y', 'py', 'd']

	def __init__(self, sigmaInitial, gaussdpp, **kwargs):
		'''
			Input:
				sigmaInitial- list; initial particle distribution 
				gaussdpp	- boll; if True, the Gaussian distribution in energy offset is used, otherwise - linear distribution
							only used with Mapclass, in the simple tracker, by default Gauss
			optional:
				order		- int; desired orer of calculation precision
				Nparticles	- int; number of particles used in the custom tracking (Not used if only MAPCLASS is used)
		'''
		
#		assert sigmaInitial != None, "sigmaInitial should be defined"

		self.sigmaInitial, self.gaussdpp, self.Nparticles = sigmaInitial, gaussdpp, kwargs.get('Nparticles', 10000)
		if 'order' in kwargs:
			self.order = kwargs.get('order')
			
		self.Map = mapclass.Map2(self.order)
						
	def update(self):
		'''
			Reads the map from "fort.18" file. Is used when the map is updated by MadX which is run externally
		'''
		self.Map = mapclass.Map2(self.order)

	def __create_set(self, index_list = [0]):
		result = [0, 0, 0, 0, 0]
		for s in index_list:
			result[s] += 1
		return result


	def __get_term(self, coordinate = 'x', set = [1, 0, 0, 0, 0]):
		if coordinate == 'x':
			return self.Map.x[set[0], set[1], set[2], set[3], set[4]]
		elif coordinate == 'px': 
			return self.Map.px[set[0], set[1], set[2], set[3], set[4]]
		elif coordinate == 'y': 
			return self.Map.y[set[0], set[1], set[2], set[3], set[4]]
		elif coordinate == 'py': 
			return self.Map.py[set[0], set[1], set[2], set[3], set[4]]
		elif coordinate == 'd': 
			return self.Map.d[set[0], set[1], set[2], set[3], set[4]]
		else:
			raise ValueError("Incorrect coordinate")

	def get_coordinate(self, coordinate, init_vector, highest_order = 2):
		'''
			Returns the coordinate-vector of the particle at the end of the lattice

			Input:
				coordinate		- string (from coordinates_core.state); coordinate to be calculated
				init_vector		- list; initial particle vector
				highest_order	- int; calculation presicion, should not be larger than 3

			**This method ignores the static component of the coordinates -> orbit 

		'''
		assert highest_order < 4

		x = None
		for s in self.state:
			if coordinate == s: x = s

		if x == None: raise ValueError

		result = 0.0

		#first order
		for i in range(5):
			result += self.__get_term(coordinate, self.__create_set([i])) * init_vector[i]
#			"""--testing--"""
#			print "R3" + str(i) + "contibution = " + self.__get_term(coordinate, self.__create_set([i])) * init_vector[i]

		if highest_order == 1: return result

		#second order
		for i in range(5):
			for j in range(5):
				result += self.__get_term(coordinate, self.__create_set([i, j])) * init_vector[i] * init_vector[j]
#					"""--testing--"""
#					print "T3" + str(i) + str(j) + "contibution = " + self.__get_term(coordinate, self.__create_set([i, j])) * init_vector[i] * init_vector[j]
		if highest_order == 2: return result

		#third order - be sure that the map is of the corresponding order
		for i in range(5):
			for j in range(5):
				for k in range(5):
					#print i, j, k
					result += self.__get_term(coordinate, self.__create_set([i, j, k])) * init_vector[i] * init_vector[j] * init_vector[k]
		if highest_order == 3: return result

		return result

	def propagate(self, Nparticles):
		'''
			Calculate the beam coordinates at the end of the lattice and stores to coordinates_core.vector_fin

			Input:
				Nparticles	- int; number of the particles in the beam
		'''

		#There is no need to store the energy deviations
		self.vector_fin = [np.zeros(self.Nparticles), np.zeros(self.Nparticles), np.zeros(self.Nparticles), np.zeros(self.Nparticles), np.zeros(self.Nparticles)]

		print("Propagating the particles..")
		for i in range(self.Nparticles):
			vector_in = [self.x_initial[i], self.px_initial[i], self.y_initial[i], self.py_initial[i], self.d_initial[i]]
			self.vector_fin[0][i] = self.get_coordinate('x', vector_in, self.order)
			self.vector_fin[1][i] = self.get_coordinate('px', vector_in, self.order)
			self.vector_fin[2][i] = self.get_coordinate('y', vector_in, self.order)
			self.vector_fin[3][i] = self.get_coordinate('py', vector_in, self.order)
			self.vector_fin[4] = self.d_initial

	def generate_init_distribution():
		self.x_initial = np.random.normal(0.0, self.sigmaInitial[0], self.Nparticles)
		self.px_initial = np.random.normal(0.0, self.sigmaInitial[1],self.Nparticles)
		self.y_initial = np.random.normal(0.0, self.sigmaInitial[2], self.Nparticles)
		self.py_initial = np.random.normal(0.0, self.sigmaInitial[3], self.Nparticles)
		self.d_initial = np.random.normal(0.0, self.sigmaInitial[4], self.Nparticles)

	def get_sigma(self, array_list, method = 0):
		'''
			Returns the sigma for a given setup in the array_list or the beam size if a single int is given

			Input:
				array_list	- list(int); indicate the sigma to be calculated giving the location of each coordinate in the coordinate_core.state
					Ex: list_of_terms = [0,1] means that we want sigma_x_x'
				method		- int (0 or 1); indicates whether Mapclass or internal tracking is used:
							0: MAPCLASS
							1: MAPCLASS + custom tracking
		'''

		if method == 1:
			if not is_local(self.x_initial):
				generate_init_distribution()
			
			self.propagate(self.Nparticles)

			if len(array_list) == 1:
				return np.sum(self.vector_fin[array_list[0]]) / self.Nparticles
			elif len(array_list) == 2:
				return np.sum(self.vector_fin[array_list[0]] * self.vector_fin[array_list[1]]) / self.Nparticles
			elif len(array_list) == 3:
				return np.sum(self.vector_fin[array_list[0]] * self.vector_fin[array_list[1]] * self.vector_fin[array_list[2]]) / self.Nparticles
			else:
				raise Exception("array_list - Order up to 3 is only supported")

		elif method == 0:
			if len(array_list) == 1:
				return self.Map.offset(self.state[array_list[0]], self.sigmaInitial, self.gaussdpp)
			elif len(array_list) == 2:
				if (array_list[0] == array_list[1]) and (array_list[0] == 2):
					return self.Map.sigma('y', self.sigmaInitial, self.gaussdpp)
				elif (array_list[0] == array_list[1]) and (array_list[0] == 0):
					return self.Map.sigma('x', self.sigmaInitial, self.gaussdpp)
				else:
					return self.Map.correlation(self.state[array_list[0]], self.state[array_list[1]], self.sigmaInitial, self.gaussdpp)
			elif len(array_list) == 3:
				return self.Map.correlation3(self.state[array_list[0]], self.state[array_list[1]], self.state[array_list[2]], self.sigmaInitial, self.gaussdpp)
			else:
				raise Exception("array_list - Order up to 3 is only supported")
		else:
			raise Exception("method should be 1 or 0")
		
		return None
	
	"""Utilities"""
	def get_obs(self, **kwargs):
		'''
			Updates the observable for the knob iteration
			Normally, beam size is used as the observable

			Input:
				observable	- type of observable
					list(int) - custom covariance checking.
				int			- custom observable
					0 - vertical beam size
					1 - horizontal beam size

			**Optionally it is possible to use any other variables, such as luminosity
		'''
		obs, method = kwargs.get('observable', 0), kwargs.get('method', 0)

		if isinstance(obs, list):
			'''calculates the custom conjugation'''
			return self.get_sigma(obs)
		elif obs == 0:
			'''default one, correspond to the vertical beam size tuning'''
			return abs(self.Map.sigma('y', self.sigmaInitial, self.gaussdpp) - self.Map.offset('y', self.sigmaInitial, self.gaussdpp) ** 2)
		elif obs == 1:
			'''Horizontal beam size tuning || only used for the waist scan'''
			return abs(self.Map.sigma('x', self.sigmaInitial, self.gaussdpp) - self.Map.offset('x', self.sigmaInitial, self.gaussdpp) ** 2)
		elif obs == 2:
			pass
			raise Exception("observable - incorrect value") 
		return None


'''Testing routine'''
if __name__ == "__main__":
	betx = 8.62817146
	bety = 0.718388011

	gamma = 2544.0313111546
	emittanceXnorm = 5.08e-6
	emittanceYnorm = 3.0e-8
	emittanceX = emittanceXnorm/gamma
	emittanceY = emittanceYnorm/gamma
	dp = 0.0008

	sigmaInitial = [np.sqrt(emittanceX * betx), np.sqrt(emittanceX / betx), np.sqrt(emittanceY * bety), np.sqrt(emittanceY / bety), dp]

	a = coordinates_core(sigmaInitial, True, order = 3, Nparticles = 10000)


