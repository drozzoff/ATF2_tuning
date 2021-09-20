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
	atf2_machine.save_knobs()

def check():

	#read the knobs from a file
	atf2_machine.knobs = map(lambda x: Knob(filename = x), ['knobs_storage/kax.knob', 'knobs_storage/kay.knob', 'knobs_storage/kex.knob'])

	#converting the knobs to the string in the suitable format for madx, and setting the calculations with them
	atf2_machine.calculation_settings['knobs_info'] = atf2_machine.save_knobs_for_madx()
	atf2_machine.calculation_settings['sext_off'] = False

	#check the knob	| range_scale to be chosen carrefully -> Mad-X error can occur
	atf2_machine.knob_check('kex', x_observables, range_scale = 1e-0)

def rescale():
	#read the knobs from a file
	atf2_machine.knobs = map(lambda x: Knob(filename = x), ['kax.knob', 'kay.knob', 'kex.knob'])

	#converting the knobs to the string in the suitable format for madx, and setting the calculations with them
	atf2_machine.calculation_settings['knobs_info'] = atf2_machine.save_knobs_for_madx()
	atf2_machine.calculation_settings['sext_off'] = False

	#rescale the knob | amplitude to be chosen carrefully (default 1.0) -> Mad-X error can occur
	atf2_machine.normalize_knob('kay', [2, 3], amplitude = 1e-3, target = 1e-7)

	#updating the MAD-X string
	atf2_machine.calculation_settings['knobs_info'] = atf2_machine.save_knobs_for_madx()

	#check the knob	| range_scale to be chosen carrefully -> Mad-X error can occur
	atf2_machine.knob_check('kay', x_observables, range_scale = 1.0)

	atf2_machine.save_knobs()

#construction()

check()

#rescale()
