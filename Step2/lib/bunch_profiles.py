# 14.10.2017: include dpp profile

import numpy as np
import orbit_mpi
import scipy.io as sio
from orbit_mpi import mpi_datatype, mpi_op
from spacecharge import Grid1D
from orbit_utils import BunchExtremaCalculator
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy.optimize import leastsq
from scipy.interpolate import UnivariateSpline

class Bunch_Profile():
	def __init__(self, grid_size_x=100, grid_size_y=100, grid_size_z=30):
		self.grid_objs = {'x': Grid1D(grid_size_x), 'y': Grid1D(grid_size_y), 'z': Grid1D(grid_size_z), 'dpp': Grid1D(grid_size_z)}
			
	def bin_bunch(self, bunch, twiss_dict):
		grid_objs = self.grid_objs
		n_part = bunch.getSize()
		sync_particle = bunch.getSyncParticle()
		beta = sync_particle.beta()
		gamma = sync_particle.gamma()
		energy = sync_particle.mass() * gamma

		# the bunch and lattice parameters are stored as attributes for later emittance calculation
		self.beta = beta
		self.gamma = gamma
		self.twiss_dict = twiss_dict
		
		dpp = np.array(map(bunch.dE, xrange(n_part))) / energy / beta**2
		coords = {}
		coords['x'] = np.array(map(bunch.x, xrange(n_part))) - twiss_dict['etax']*dpp
		coords['y'] = np.array(map(bunch.y, xrange(n_part))) - twiss_dict['etay']*dpp
		coords['z'] = np.array(map(bunch.z, xrange(n_part)))
		coords['dpp'] = dpp

		comm = bunch.getMPIComm()		
		self.grid_arrs_dict = {}
		self.grid_arrs_norm_dict = {}
		for u in ['x', 'y', 'z', 'dpp']:
			self._bin_coordinate(grid_objs[u], coords[u], comm)
			self.grid_arrs_dict[u] = self._get_grid_arrs(grid_objs[u])
			self.grid_arrs_norm_dict[u] = self._get_grid_normalized_amplitude(self.grid_arrs_dict[u])		
		return self.grid_arrs_dict
	
	def save_bunchprofile_to_matfile(self, matfile_out):
		rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
		if not rank:
			with open(matfile_out,'w') as f:
				sio.savemat(f, self.grid_arrs_dict)

	def _bin_coordinate(self, grid, u, comm):
		grid.setZero()
		Min = orbit_mpi.MPI_Allreduce(min(u),mpi_datatype.MPI_DOUBLE,mpi_op.MPI_MIN,comm)
		Max = orbit_mpi.MPI_Allreduce(max(u),mpi_datatype.MPI_DOUBLE,mpi_op.MPI_MAX,comm)
		grid_size = grid.getSizeZ()
		delta = (Max-Min)/grid_size
		Min -= delta*1.5
		Max += delta*1.5
		grid.setGridZ(Min, Max) 
		map(lambda i: grid.binValue(1, u[i]), xrange(len(u)))
		grid.synchronizeMPI()

	def _Gauss(self,x,x0,a,sigma):
	    return a/(np.sqrt(2*np.pi)*np.abs(sigma))*exp(-(x-x0)**2/(2*sigma**2))

	def _GaussianFit(self, x, y):
		mean = sum(x*y)/sum(y)
		sigma = np.sqrt(sum(y*(x-mean)**2)/sum(y))
		dx = x[1]-x[0]
		amplitude = sum(y)*dx
		popt,pcov = curve_fit(self._Gauss,x,y,p0=[mean,amplitude,sigma])
		return popt

	def transverse_emittances_Gauss(self):
		amplitude_normalized = {}
		epsn = {}
		for u in ['x', 'y']:
			popt = self._GaussianFit(*self.grid_arrs_norm_dict[u])
			amplitude_normalized[u] = popt[1]
			epsn[u] = popt[-1]**2 * self.gamma*self.beta/self.twiss_dict['beta'+u]
		return epsn, amplitude_normalized

	def dpp_Gauss_fitted(self):
		popt = self._GaussianFit(*self.grid_arrs_norm_dict['dpp'])
		return popt[2]

	def dpp_from_FWHM(self):
		x, y = self.grid_arrs_norm_dict['dpp']
		spline = UnivariateSpline(x, y-np.max(y)/2, s=0)
		r1, r2 = spline.roots() 
		return (r2-r1)/2.3548

	def _get_grid_arrs(self, grid):
		x = np.linspace(grid.getMinZ(), grid.getMaxZ(), grid.getSizeZ()+1)
		y = map(grid.getValue, x)
		return x, y

	def _get_grid_normalized_amplitude(self, (x, y)):
		# this works only for uniform x ...
		dx = x[1]-x[0]
		y /= np.sum(y)*dx
		return x, y

	def transverse_emittances_combinedDoubleGauss(self):
		def _doubleGauss(x,x0,a1,sigma1,a2,sigma2):
			return abs(a1)/(np.sqrt(2*np.pi)*np.abs(sigma1))*np.exp(-(x-x0)**2/(2*sigma1**2)) + \
				   abs(a2)/(np.sqrt(2*np.pi)*np.abs(sigma2))*np.exp(-(x-x0)**2/(2*sigma2**2))

		def err(p, x, y):
			return _doubleGauss(x, *p) - y

		def err_global(p, x_pos, y_pos, x_ampl, y_ampl):
			# p is now: a1, x0, y0, sigma1x, sigma1y, sigma2x, sigma2y
			px = p[1], abs(p[0]), p[3], 1-abs(p[0]), p[5]
			py = p[2], abs(p[0]), p[4], 1-abs(p[0]), p[6]
			errx = err(px, x_pos, x_ampl)
			erry = err(py, y_pos, y_ampl)
			return np.concatenate((errx, erry))

		def _combinedDoubleGaussianFit( (x_pos, x_ampl), (y_pos, y_ampl) ):
			px = self._GaussianFit(x_pos,x_ampl)
			mean_x = px[0]
			sigma_x = px[2]
			py = self._GaussianFit(y_pos,y_ampl)
			mean_y = py[0]
			sigma_y = py[2]
			p_global = [0.9, mean_x, mean_y, sigma_x, sigma_y, 2*sigma_x, 2*sigma_y]
			popt, ier = leastsq(err_global, p_global, args=(x_pos, y_pos, x_ampl, y_ampl), ftol=1e-12)
			return popt

		popt = _combinedDoubleGaussianFit(self.grid_arrs_norm_dict['x'], self.grid_arrs_norm_dict['y'])
		epsn_1 = {'x': popt[3]**2 * self.gamma*self.beta/self.twiss_dict['betax'],
				  'y': popt[4]**2 * self.gamma*self.beta/self.twiss_dict['betay']}
		epsn_2 = {'x': popt[5]**2 * self.gamma*self.beta/self.twiss_dict['betax'],
				  'y': popt[6]**2 * self.gamma*self.beta/self.twiss_dict['betay']}
		ampl_1 = abs(popt[0])
		ampl_2 = 1-abs(popt[0])
		return epsn_1, epsn_2, ampl_1, ampl_2



