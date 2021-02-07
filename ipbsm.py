#!/usr/bin/env python2.7

import numpy as np

from numpy.random import uniform

class IPBSM:
	"""
		This class is used to define the properties of the IP beam size monitor at ATF2

		Typical use is to incude the errors to the measurement depeneding on the current regime
	"""
	def __init__(self):
		self.properties = {
			"wire":	{"theta": None, "error": 8e-7},
			"2":	{"theta": 0.03490658503988659, "error": 1e-7},
			"6.4":	{"theta": 0.1117010721276371, "error": 1e-7},
			"8":	{"theta": 0.13962634015954636, "error": 1e-7},
			"30":	{"theta": 0.5235987755982988, "error": 2e-8},
			"174":	{"theta": 3.036872898470133, "error": 8e-9}
		}

		pass

	def __str__(self):
		print self.__dict__

	__repr__ = __str__

	def add_error(self, input_size, mode = None):
		'''
			** Should bre rewritten as a decorator
			Adds the error to the input beam size based on the absolute value and returnes the beam size with error included.

			
			Input:
				input_size	- double; input beam size, is assumed to be squared
			
			Output:
				beam size with an error, bot_error value, top_error value
			Note:
				Errors distributed according to the uniform distribution

				Typical error values of the IPBSM are used:

				Wire scanner:
					1.4+ micron : 800 nm

				8 degree mode:
					360 nm - 1.4 micron : 100 nm

				30 degree mode:
					100 nm - 360 nm : 20 nm

				174 degree mode:
					20 nm - 100 nm : 8 nm

				***
					Apart from the applying an error, a bottom level of the beam size of 20 nm is applied
				***

				mode: 
					None	- error is applied based on the actual beam size
					"2"		- 2 degree mode error [No data]
					"6.4"	- 6.4 degree mode error [No data]
					"8"		- 8 degree mode error
					"30"	- 30 degree mode error
					"174"	- 174 degree mode error
				
				***
		''' 

		beam_size = np.sqrt(input_size)
		if mode == None:
			if beam_size > 1.4e-6:
				return (beam_size + np.random.uniform(-8e-7, 8e-7))**2#, 8e-7

			if beam_size < 1.4e-6 and beam_size > 3.6e-7:
				return (beam_size + np.random.uniform(-1e-7, 1e-7))**2#, 1e-7

			if beam_size < 3.6e-7 and beam_size > 1e-7:
				return (beam_size + np.random.uniform(-2e-8, 2e-8))**2#, 2e-8

			if beam_size < 1e-7 and beam_size > 2e-8:
				return (beam_size + np.random.uniform(-8e-9, 8e-9))**2#, 8e-9

			if beam_size < 2e-8:
				return (2e-8 + np.random.uniform(0, 8e-9))**2
		
		else: 
			try:
				err = self.properties[mode]['error']
				return (beam_size + uniform(-err, err))**2#, err
			except:
				raise AttributeError(mode + " - incorrect mode")

	def transform_to_modulation(self, input_size, mode):
		'''
			Transforms the input beam size to the modulation (for a given IPBSM degree mode) and returns its value

			Input:
				input_size	- float; input beam size
				mode		- string; measurement mode used
		'''
		
#		if self.__check_mode(mode) == False:
#			raise ValueError("mode - non existing value")
		
		beam_size = np.sqrt(input_size)

		if mode == "wire":
			raise Exception("Not possible to use modulation with wire scanners")

		theta = self.properties[mode]["theta"]

		Modulation = np.abs(np.cos(theta)) * np.exp(-2 * (2 * np.pi / 532e-9 * (np.sin(theta / 2)) * beam_size )**2)

		return Modulation

#		***Modulation reduction factors***
#			to be written...

		pass


if __name__ == "__main__":
	a = IPBSM()
