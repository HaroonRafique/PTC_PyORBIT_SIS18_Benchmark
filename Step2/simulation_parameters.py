import numpy as np

# Switch for SIS18 parameters
Step1to6 = 0

# Constants
beta = 0.15448
gamma = 1.012149995
m = 1.2
TransverseCut = 5
circumference= 216.7199935

# Beam Parameters
if Step1to6:
	intensity=2.95e10	# DeltaQ = 0.1  Z=40.3m t=3472.7ns
	blength_rms = (beta*299792458*3472.7e-9)/4.	# = 40.206868 m
	sig_z = (3472.7e-9)/4.
	RestoringForce = -1.951E-11
else:
	intensity=1.95e+9	# DeltaQ = 0.1  Z=2.69m t=231.51ns
	blength_rms = (beta*299792458*231.51e-9)/4. # = 2.680419244 m
	sig_z = (231.51e-9)/4.
	RestoringForce = -4.38975e-09

# Assume Geometric epsn_x = (beta*gamma)*epsn_g_x
epsn_x=(beta*gamma)*(12.57e-6)/4			# beta*gamma*e_g/4 = 4.91E-7
epsn_y=(beta*gamma)*(9.30e-6)/4				# beta*gamma*e_g/4 = 3.635E-7
dpp_rms = 2.5e-4/3.

# Simulation Parameters
n_macroparticles = int(100)
macrosize = intensity/float(n_macroparticles)
turns_max = 1024
turns_update = range(-1, turns_max, 100)
turns_print = range(-1, turns_max, 1)
rf_voltage=0.0

parameters = {
	'LongitudinalJohoParameter': m,
	'TransverseCut': TransverseCut,
	'LongitudinalCut': 2.4,
	'circumference':circumference,
	'intensity': intensity,
	'blength_rms': blength_rms,
	'sig_z':sig_z,
	'epsn_x': epsn_x,
	'epsn_y': epsn_y,
	'dpp_rms': dpp_rms,
	'n_macroparticles': n_macroparticles,
	'macrosize': macrosize,
	'turns_max': turns_max,
	'turns_update': turns_update,
	'turns_print': turns_print,
	'rf_voltage': rf_voltage
}

switches = {
	'Horizontal': 1,
	'SliceBySlice': 0,
	'Frozen': 1,
	'MinPathLength': 1E-8,
	'RestoringForce': RestoringForce,
	'InitialParticleTransversePosition':5E-3,	# Used for single particle sims
	'InitialDistnSigma':4						# Used for bunch sims
}
