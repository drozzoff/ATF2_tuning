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

#generate the Response Matrix
atf2_machine.build_response_matrix(1e-6, variables = x_shifts, relative = False, sext_off = False, list_of_terms = x_observables)

#Extract the knobs from the precalculated Response Matrix
knobs = atf2_machine.construct_knobs(x_shifts_knobs_names, 'hor shift knobs')

for x in knobs:
	print x