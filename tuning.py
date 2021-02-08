#!/usr/bin/env python2.7

import numpy as np

from atf2tuning import ATF2Tuning, atf2_generate_errors
from machine import AbstractMachine, Knob


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

		plt.show()

betx = 8.62817146
bety = 0.718388011

gamma = 2544.0313111546
emittanceXnorm = 5.08e-6
emittanceYnorm = 3.0e-8
emittanceX = emittanceXnorm/gamma
emittanceY = emittanceYnorm/gamma
dp = 0.0008

sigmaInitial = [np.sqrt(emittanceX * betx), np.sqrt(emittanceX / betx), np.sqrt(emittanceY * bety), np.sqrt(emittanceY / bety), dp]
gaussdpp = True

atf2_machine = AbstractMachine(sigmaInitial, gaussdpp, tuning_order = 5, Nparticles = 100000, method = 0, name = "ATF2")

knobs_loc, knobs_list = "knobs_storage/", ["kax.knob", "kay.knob", "kex.knob", "coup_yx.knob", "coup_ypx.knob", "key.knob"]
atf2_machine.knobs = list(map(lambda x: Knob(filename = knobs_loc + x), knobs_list))

atf2_generate_errors()
atf2 = ATF2Tuning(atf2_machine, ipbsm_errors = False, initial_errors = "initial_errors.madx")

#tuning the linear knobs
hor_knobs = ["kax", "kex"]
vert_knobs = ["kay", "key", "coup_ypx"]

iteration_list = vert_knobs
#atf2.read_setup()
#atf2.measure_orbit(plot = plot_orbit)
atf2.correct_orbit(plot = plot_orbit)

#waist scan
atf2.tune("kqf1ff", 1, knob_range = [1.661885988 * 0.99, 1.661885988 * 1.01], additive = False, plot = iteration_plot)

atf2.tune("kqd0ff", knob_range = [(-2.857321177) * 1.01, (-2.857321177) * 0.99], additive = False, plot = iteration_plot)
__ = map(lambda x: atf2.tune(x, sext_off = False), iteration_list)	#traditional tuning with linear knobs
atf2.save_log()
atf2.save_setup()


atf2.tune("QD0TT_knob", knob_range = [-1.0, 1.0], additive = False, sext_off = False)
atf2.tune("koct2", knob_range = [-7309.9998, 7309.9998], sext_off = False, oct_off = False)

#	print atf2.knobs_preset
