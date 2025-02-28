import math
import sys
from itertools import chain
import numpy as np
import csv
import random
import orbit_mpi

from bunch import Bunch
from orbit.injection.joho import JohoLongitudinal
from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import WaterBagDist2D, GaussDist2D, KVDist2D
from orbit.utils.consts import mass_proton, speed_of_light, pi


class JohoNormalized2D:

	def __init__(self, order):
		"""
		Constructor.
		"""
		self.order = order
		self.orderinv = 1. / order
		self.emitlim = (1. + order) / 0.5
		self.__initialize()
		
	def __initialize(self):
		self.emit = self.emitlim * 2./(1. + self.order)
		self.poslength = math.sqrt(self.emitlim)
		self.momlength = math.sqrt(self.emitlim)
		self.emitrms = 0.5 * self.emitlim/(1. + self.order)

	def getCoordinates(self):
		s1 = random.random()
		s2 = random.random()
		a = math.sqrt(1 - pow(s1, self.orderinv))
		al = 2. * math.pi * s2
		u = a * math.cos(al)
		v = a * math.sin(al)
		pos = self.poslength * u
		mom = self.momlength * v		
		return (pos,mom)

class LongitudinalJohoDistributionSingleHarmonic():

	def __init__(self, Parameters_dict, Joho_order):
		"""
		Constructor.
		"""	
		self.beta = Parameters_dict['beta']
		circumference = Parameters_dict['circumference']
		self.harmonic_number = Parameters_dict['harmonic_number']
		gamma_transition = Parameters_dict['gamma_transition']
		gamma = Parameters_dict['gamma']
		bunch_length = Parameters_dict['bunch_length']
		self.phi_s = Parameters_dict['phi_s']
		self.rf_voltage = Parameters_dict['rf_voltage']
		self.LongitudinalCut = Parameters_dict['LongitudinalCut']
	
		self.rf_frequency = speed_of_light * self.beta / circumference * self.harmonic_number
		self.slip_factor = (1 / gamma_transition**2 - 1 / gamma**2)	
		self.energy = Parameters_dict['energy']
		self.sigma_dphi = self.rf_frequency * (bunch_length/2) * pi
		self.sigma_dE = np.sqrt(self.beta**2 * self.energy / self.slip_factor / self.harmonic_number * self.rf_voltage / pi * (np.cos(self.sigma_dphi + self.phi_s) - np.cos(self.phi_s) + self.sigma_dphi * np.sin(self.phi_s)))
		
		self.H_max = (self.rf_voltage/(2*pi) * (np.cos(self.LongitudinalCut * self.sigma_dphi + self.phi_s) - np.cos(self.phi_s) + (self.LongitudinalCut * self.sigma_dphi)*np.sin(self.phi_s)))
		self.distribution = JohoNormalized2D(Joho_order)
		
	def is_inside_limiting_countour(self, phi, dE):
		H = abs(self.harmonic_number * self.slip_factor / (2*self.beta**2 * self.energy) * dE**2 + self.rf_voltage/(2*pi)*(np.cos(phi) - np.cos(self.phi_s) + (phi-self.phi_s)*np.sin(self.phi_s))) 
		return H <= abs(self.H_max)

	def getCoordinates(self):
		while True:
			phi, dE = self.distribution.getCoordinates()
			phi *= self.sigma_dphi
			dE *= self.sigma_dE
			if self.is_inside_limiting_countour(phi, dE): 
				break
		return phi, dE
		

def generate_initial_distribution(parameters, output_file = 'Input/ParticleDistribution.in', summary_file = 'Input/ParticleDistribution_summary.txt', outputFormat='Orbit'):

	# twiss containers
	twissX = TwissContainer(alpha = parameters['alphax0'], beta = parameters['betax0'], emittance = parameters['epsn_x'] / parameters['gamma'] / parameters['beta'])
	twissY = TwissContainer(alpha = parameters['alphay0'], beta = parameters['betay0'], emittance = parameters['epsn_y'] / parameters['gamma'] / parameters['beta'])
	dispersionx = {'etax0': parameters['etax0'], 'etapx0': parameters['etapx0']}
	dispersiony = {'etay0': parameters['etay0'], 'etapy0': parameters['etapy0']}
	closedOrbitx = {'x0': parameters['x0'], 'xp0': parameters['xp0']} 
	closedOrbity = {'y0': parameters['y0'], 'yp0': parameters['yp0']} 

	# initialize particle arrays
	x = np.zeros(parameters['n_macroparticles'])
	xp = np.zeros(parameters['n_macroparticles'])
	y = np.zeros(parameters['n_macroparticles'])
	yp = np.zeros(parameters['n_macroparticles'])
	phi = np.zeros(parameters['n_macroparticles'])
	dE = np.zeros(parameters['n_macroparticles'])

	# building the distributions
	Transverse_distribution = GaussDist2D(twissX, twissY, cut_off=parameters['TransverseCut'])
	Longitudinal_distribution = LongitudinalJohoDistributionSingleHarmonic(parameters, parameters['LongitudinalJohoParameter'])

	if orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD) == 0:
		fid = open(output_file,"w")
		csv_writer = csv.writer(fid, delimiter=' ')
		for i in range(parameters['n_macroparticles']):
			(phi[i], dE[i]) = Longitudinal_distribution.getCoordinates()
			(x[i], xp[i], y[i], yp[i]) = Transverse_distribution.getCoordinates()
			x[i] += closedOrbitx['x0']
			xp[i] += closedOrbitx['xp0']
			y[i] += closedOrbity['y0']
			yp[i] += closedOrbity['yp0']
			dpp = dE[i] / (parameters['energy']) / parameters['beta']**2
			x[i] += dpp * dispersionx['etax0']
			xp[i] += dpp * dispersionx['etapx0']	
			y[i] += dpp * dispersiony['etay0']
			yp[i] += dpp * dispersiony['etapy0']	
		
			if outputFormat == 'Orbit':
				x[i] *= 1000.
				xp[i] *= 1000.
				y[i] *= 1000.
				yp[i] *= 1000.
				dE[i] /= 1.e9		
				csv_writer.writerow([x[i], xp[i], y[i], yp[i], phi[i], dE[i]])
		#	else:
				# still need to convert from phi to z!!
				#csv_writer.writerow([x[i], xp[i], y[i], yp[i], z[i], dE[i]])		
		fid.close()

		fid = open(summary_file, 'w')
		parameter_list = ['circumference', 'rf_voltage', 'phi_s', 'harmonic_number', 'gamma_transition', 'n_macroparticles', 'energy', 'gamma', 'bunch_length', 'LongitudinalCut', 'LongitudinalJohoParameter', 'x0', 'xp0', 'betax0', 'alphax0', 'etax0', 'etapx0', 'y0', 'yp0', 'betay0', 'alphay0', 'etay0', 'etapy0', 'epsn_x', 'epsn_y', 'TransverseCut']
		for key in parameter_list:
			fid.write(key + ' = ' + str(parameters[key]) + '\n')
		fid.close()

		print '\nCreated particle distribution with ' + str(parameters['n_macroparticles']) + ' macroparticles into file: ', output_file
	
	return output_file


