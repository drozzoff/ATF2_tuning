#!/usr/bin/env python2.7

import numpy as np

from machine import AbstractMachine, Knob

#initial beam parameters
betx = 8.62817146
bety = 0.718388011

gamma = 2544.0313111546
emittanceXnorm = 5.08e-6
emittanceYnorm = 3.0e-8
emittanceX = emittanceXnorm/gamma
emittanceY = emittanceYnorm/gamma
dp = 0.0008

#	These variables should exist in MadX scripts
#Linear knobs 
y_shifts = ["SF6DY", "SF5DY", "SD4DY", "SF1DY", "SD0DY"]
y_shifts_knobs_names = [ "coup_yx", "coup_ypx", "key"]
y_shifts_knobs_names_full = [ "coup_yx", "coup_ypx", "key", "coup_pypx"]
y_observables = [[2, 0], [2, 1], [2, 4]]
y_observables_full = [[2, 0], [2, 1], [2, 4], [3, 1]]

x_shifts = ["SF6DX", "SF5DX", "SD4DX", "SF1DX", "SD0DX"]
x_shifts_knobs_names = ["kax", "kay", "kex"]
x_shifts_knobs_names_full = ["kbetx", "kbety", "kax", "kay", "kex"]
x_observables = [[0, 1], [2, 3], [0, 4]]
x_observables_full = [[0, 0], [2, 2], [0, 1], [2, 3], [0, 4]]
#Nonlinear knobs
kn_iteration = ['KSF6KN', 'KSF5KN', 'KSD4KN', 'KSF1KN', 'KSD0KN']
kn_knobs_names = ["k324", "k346"]
kn_observables = [[2, 1, 3], [2, 3, 4]]

ks_iteration = ['KSK4KS', 'KSK3KS', 'KSK2KS', 'KSK1KS']
ks_knobs_names = ["k322", "k326", "k344", "k366"]
ks_observables = [[2, 1, 1], [2, 1, 4], [2, 3, 3], [2, 4, 4]]


sigmaInitial = [np.sqrt(emittanceX * betx), np.sqrt(emittanceX / betx), np.sqrt(emittanceY * bety), np.sqrt(emittanceY / bety), dp]
gaussdpp = True

#initialize the beamline
atf2_machine = AbstractMachine(sigmaInitial, gaussdpp, order = 2, Nparticles = 100000, method = 0, name = "ATF2")


def construction():

	#generate the Response Matrix
	atf2_machine.build_response_matrix(1e-6, variables = x_shifts, relative = False, list_of_terms = x_observables, sext_off = False)

	#Extract the knobs from the precalculated Response Matrix
	knobs = atf2_machine.construct_knobs(x_shifts_knobs_names, 'hor shift knobs')

	#saving knobs
#	atf2_machine.save_knobs()

def check():

	#read the knobs from a file
	atf2_machine.knobs = map(lambda x: Knob(filename = x), ['knobs_storage/kax.knob', 'knobs_storage/kay.knob', 'knobs_storage/kex.knob'])

	#converting the knobs to the string in the suitable format for madx, and setting the calculations with them
	atf2_machine.calculation_settings['knobs_info'] = atf2_machine.save_knobs_for_madx()
	atf2_machine.calculation_settings['sext_off'] = False

	#check the knob	| range_scale to be chosen carrefully -> Mad-X error can occur
	iter_data = atf2_machine._knob_check('kay', x_observables)
	print iter_data

def rescale():
	#read the knobs from a file
	atf2_machine.knobs = map(lambda x: Knob(filename = x), ['knobs_storage/kax.knob', 'knobs_storage/kay.knob', 'knobs_storage/kex.knob'])

	#converting the knobs to the string in the suitable format for madx, and setting the calculations with them
	atf2_machine.calculation_settings['knobs_info'] = atf2_machine.save_knobs_for_madx()
	atf2_machine.calculation_settings['sext_off'] = False
	
	#check the knob	| range_scale to be chosen carrefully -> Mad-X error can occur
	atf2_machine.knob_check('kay', x_observables, range_scale = 1.0)

	#rescale the knob | amplitude to be chosen carrefully (default 1.0) -> Mad-X error can occur
	atf2_machine.normalize_knob('kay', [2, 3], amplitude = 1e2, target = 1e-7)

	#updating the MAD-X string
	atf2_machine.calculation_settings['knobs_info'] = atf2_machine.save_knobs_for_madx()

	#check the knob	| range_scale to be chosen carrefully -> Mad-X error can occur
	atf2_machine.knob_check('kay', x_observables, range_scale = 1.0)

def simulate_tuning():
	from atf2tuning import ATF2Tuning, atf2_generate_errors


	#plotting routines
	def iteration_plot(knob_range, observables, fit_function = None, **kwargs):
		'''plotting routine'''
		import matplotlib.pyplot as plt
		x = np.linspace(knob_range[0], knob_range[-1], num = 100)
		
		with plt.style.context(['science', 'ieee']):
			plt.ylabel(kwargs.get('ylabel', "M"))
			plt.xlabel(kwargs.get('xlabel', "Knob Amplitude"))

			diff = (max(knob_range) - min(knob_range)) * 0.025
			plt.xlim(min(knob_range) - diff, max(knob_range) + diff)
			plt.plot(knob_range, observables, 'o', label = "measurement", color = "black", markersize = 3)
			if fit_function != None:
				fit_func = list(map(lambda a: fit_function(a), x))
				plt.ylim(0.0, max([max(fit_func), max(observables)]) * 1.025)
				plt.plot(x, fit_func, label = "Fit", color = "red")
			else:
				plt.ylim(0.0, max(observables) * 1.025)
			plt.legend()

			if 'filename' in kwargs:
				plt.savefig(kwargs['filename'])
			plt.show()

	def plot_orbit(bpm_list, plane = 'x', **kwargs):
		'''bpm_list is dict'''
		import matplotlib.pyplot as plt
		assert (plane == 'x') or (plane == 'y'), "Incorrect plane"

		bpm_list_cleared = {}
		for bpm in bpm_list:
			if bpm_list[bpm]['status'] == "on": bpm_list_cleared[bpm] = bpm_list[bpm]

		with plt.style.context(['science', 'ieee']):
			plt.xlabel('s [m]')
			plt.ylabel(kwargs.get('ylabel', "x"))
			key_order = sorted(bpm_list_cleared, key = lambda x: bpm_list_cleared[x]['s'])

			plt.plot(list(map(lambda x: bpm_list_cleared[x]['s'], key_order)), list(map(lambda x: bpm_list_cleared[x][plane], key_order)), '-o', markersize = 3)


	#read the knobs from a file
	knobs_loc, knobs_list = "knobs_storage/", ["kax.knob", "kay.knob", "kex.knob", "coup_yx.knob", "coup_ypx.knob", "key.knob"]
	atf2_machine.knobs = map(lambda x: Knob(filename = knobs_loc + x), knobs_list)

	#generate rrors for MAD-X
	atf2_generate_errors()
	atf2 = ATF2Tuning(atf2_machine, ipbsm_errors = False, initial_errors = "initial_errors.madx", mode = "wire")

	#Correct orbit with MAD internal routine
	atf2.correct_orbit(plot = plot_orbit)
	return
	#waist scan
	atf2.tune("kqf1ff", 1, knob_range = [1.661885988 * 0.99, 1.661885988 * 1.01], additive = False)
	atf2.tune("kqd0ff", knob_range = [(-2.857321177) * 1.01, (-2.857321177) * 0.99], additive = False)

	#tuning with the predesigned tuning knobs
	__ = map(lambda x: atf2.tune(x, sext_off = False), ["kay", "key", "coup_ypx"])	#traditional tuning with linear knobs
	atf2.save_log()
	atf2.save_setup()

def octupole_alignment():
	from atf2tuning import ATF2Tuning, atf2_generate_errors
	from chainmap import ChainMap

	def iteration_plot(knob_range, observables, fit_function = None, **kwargs):
		'''plotting routine'''
		import matplotlib.pyplot as plt
		x = np.linspace(knob_range[0], knob_range[-1], num = 100)
		print knob_range, observables
		with plt.style.context(['science', 'ieee']):
			plt.ylabel(kwargs.get('ylabel', "M"))
			plt.xlabel(kwargs.get('xlabel', "Knob Amplitude"))

			diff = (max(knob_range) - min(knob_range)) * 0.025
			plt.xlim(min(knob_range) - diff, max(knob_range) + diff)
			plt.plot(knob_range, observables, 'o', label = "measurement", color = "black", markersize = 3)
			y_left, y_right, x_right, x_left = None, None, None, None

			if 'ylim' in kwargs:
				y_left, y_right = kwargs.get('ylim')
			if 'xlim' in kwargs:
				x_left, x_right = kwargs.get('xlim')

			if fit_function != None:
				fit_func = map(lambda a: fit_function(a), x)
				if 'ylim' in kwargs:
					plt.ylim(y_left, y_right)
				else:
					diff_y = (max([max(fit_func), max(observables)]) - min([min(fit_func), min(observables)])) * 0.025
					plt.ylim(min([min(fit_func), min(observables)]) - diff_y, max([max(fit_func), max(observables)]) + diff_y)
				plt.plot(x, fit_func, label = "Fit", color = "red")
			else:
				if 'ylim' in kwargs:
					plt.ylim(y_left, y_right)
				else:
					diff_y = (max(observables) - min(observables)) * 0.025
					plt.ylim(min(observables) - diff_y, max(observables) + diff_y)
			plt.legend()

			if 'filename' in kwargs:
				plt.savefig(kwargs['filename'])
			plt.show()

#	def oct_plot(*args, **kwargs):
#		iteration_plot(*args, **kwargs, xlabel = "")

	atf2_machine.knobs = [Knob(filename = 'knobs_storage/kay.knob')]

	atf2 = ATF2Tuning(atf2_machine, ipbsm_errors = False, mode = "174")

	atf2.oct2_alignment(offset_points = 7, knob_range = [-1.0, 1.0], bba_plot = lambda *args, **kwargs: iteration_plot(*args, **ChainMap({'xlabel': "OCT2 offset", 'ylabel': "Waist shift", 'filename': "test.pdf"}, kwargs)))

#construction()

#check()

#rescale()

#simulate_tuning()

octupole_alignment()