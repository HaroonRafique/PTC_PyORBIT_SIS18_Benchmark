#----------------------------------------------
# Simulation Switches
#----------------------------------------------
from simulation_parameters import switches as s

slicebyslice = s['SliceBySlice']        # 2.5D space charge
frozen = s['Frozen']                    # Frozen space charge
if frozen:
        slicebyslice=0

horizontal = s['Horizontal']            # Horizontal Poincare Distn

import math
import sys
import time
import orbit_mpi
import timeit
import numpy as np
import scipy.io as sio
import os

# utils
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from orbit.utils.consts import mass_proton, speed_of_light, pi

# bunch
from bunch import Bunch
from bunch import BunchTwissAnalysis, BunchTuneAnalysis
from orbit.bunch_utils import ParticleIdNumber

# diagnostics
from orbit.diagnostics import TeapotStatLatsNode, TeapotMomentsNode, TeapotTuneAnalysisNode
from orbit.diagnostics import addTeapotDiagnosticsNodeAsChild
from orbit.diagnostics import addTeapotMomentsNodeSet, addTeapotStatLatsNodeSet

# PTC lattice
from libptc_orbit import *
from ext.ptc_orbit import PTC_Lattice
from ext.ptc_orbit import PTC_Node
from ext.ptc_orbit.ptc_orbit import setBunchParamsPTC, readAccelTablePTC, readScriptPTC
from ext.ptc_orbit.ptc_orbit import updateParamsPTC, synchronousSetPTC, synchronousAfterPTC
from ext.ptc_orbit.ptc_orbit import trackBunchThroughLatticePTC, trackBunchInRangePTC
from orbit.aperture import TeapotApertureNode

# transverse space charge
if frozen:
        from orbit.space_charge.analytical import scAccNodes, scLatticeModifications
        from spacecharge import SpaceChargeCalcAnalyticGaussian
        from spacecharge import GaussianLineDensityProfile
if slicebyslice:
        #PIC
        from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications
        from spacecharge import SpaceChargeCalcAnalyticGaussian
        from spacecharge import InterpolatedLineDensityProfile


from lib.output_dictionary import *
from lib.particle_output_dictionary import *
from lib.pyOrbit_GenerateInitialDistribution2 import *
from lib.save_bunch_as_matfile import *
from lib.pyOrbit_LinearRestoringForce import *

print "Start ..."
comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
rank = orbit_mpi.MPI_Comm_rank(comm)

#----------------------------------------------
# Create folder structure
#----------------------------------------------
from lib.mpi_helpers import mpi_mkdir_p
mpi_mkdir_p('input')
mpi_mkdir_p('output')
mpi_mkdir_p('lost')

#----------------------------------------------
# Generate Lattice (MADX + PTC)
#----------------------------------------------
if not rank:
	os.system("/afs/cern.ch/eng/sl/MAD-X/pro/releases/5.02.00/madx-linux64 < Input/SIS18.madx")
orbit_mpi.MPI_Barrier(comm)

#----------------------------------------------
# Initialize a Teapot-Style PTC lattice
#----------------------------------------------
PTC_File = "SIS_18_BENCHMARK.flt"
Lattice = PTC_Lattice("MACHINE")
Lattice.readPTC(PTC_File)
readScriptPTC('Input/time.ptc')

paramsDict = {}
paramsDict["length"]=Lattice.getLength()/Lattice.nHarm

#----------------------------------------------
# Add apertures
#----------------------------------------------
position = 0
for node in Lattice.getNodes():
	myaperturenode = TeapotApertureNode(1, 10, 10, position)
	node.addChildNode(myaperturenode, node.ENTRANCE)
	node.addChildNode(myaperturenode, node.BODY)
	node.addChildNode(myaperturenode, node.EXIT)
	position += node.getLength()

#----------------------------------------------
# Add the main bunch and lost particles bunch
#----------------------------------------------
bunch = Bunch()
setBunchParamsPTC(bunch)

from simulation_parameters import parameters as p
p['harmonic_number'] = Lattice.nHarm 
p['phi_s']           = 0
p['gamma']           = bunch.getSyncParticle().gamma()
p['beta']            = bunch.getSyncParticle().beta()
p['energy']          = 1e9 * bunch.mass() * bunch.getSyncParticle().gamma()
p['bunch_length'] = p['blength_rms']/speed_of_light/bunch.getSyncParticle().beta()*4
kin_Energy = bunch.getSyncParticle().kinEnergy()
print 'Energy of particle = ', p['energy']
print 'Kinetic Energy of particle = ', kin_Energy

if horizontal: 
	Particle_distribution_file = generate_initial_5mm_distributionH(s['InitialParticleTransversePosition'], 0, p, Lattice, output_file='input/ParticleDistribution.in', summary_file='input/ParticleDistribution_summary.txt')
else:
	Particle_distribution_file = generate_initial_5mm_distributionV(s['InitialParticleTransversePosition'], 0, p, Lattice, output_file='input/ParticleDistribution.in', summary_file='input/ParticleDistribution_summary.txt')

bunch_orbit_to_pyorbit(paramsDict["length"], kin_Energy, Particle_distribution_file, bunch, p['n_macroparticles'] + 1) #read in only first N_mp particles.
bunch.addPartAttr("macrosize")
map(lambda i: bunch.partAttrValue("macrosize", i, 0, p['macrosize']), range(bunch.getSize()))
ParticleIdNumber().addParticleIdNumbers(bunch) # Give them unique number IDs
bunch.dumpBunch("input/mainbunch_start.dat")
saveBunchAsMatfile(bunch, "output/mainbunch_-000001")

lostbunch = Bunch()
bunch.copyEmptyBunchTo(lostbunch)
lostbunch.addPartAttr('ParticlePhaseAttributes')
lostbunch.addPartAttr("LostParticleAttributes")
paramsDict["lostbunch"]=lostbunch
paramsDict["bunch"]= bunch

#----------------------------------------------------
# Add space charge nodes - FROZEN
#----------------------------------------------------
if frozen:
        print '\nSetting up the space charge calculations ...'
        # Make a SC solver using frozen potential
        # ~ sc_path_length_min = 1E-8
        sc_path_length_min = s['MinPathLength']
        LineDensity=GaussianLineDensityProfile(p['blength_rms'])
        sc_params1 = {'intensity': p['intensity'], 'epsn_x': p['epsn_x'], 'epsn_y': p['epsn_y'], 'dpp_rms': p['dpp_rms'], 'LineDensity': LineDensity}
        space_charge_solver1 = SpaceChargeCalcAnalyticGaussian(*[sc_params1[k] for k in ['intensity','epsn_x','epsn_y','dpp_rms','LineDensity']])
        print dir(scLatticeModifications)
        sc_nodes1 = scLatticeModifications.setSCanalyticalAccNodes(Lattice, sc_path_length_min, space_charge_solver1)
        print 'Installed %i space charge nodes'%(len(sc_nodes1))

#----------------------------------------------------
# Add space charge nodes - SliceBySlice
#----------------------------------------------------
if slicebyslice:
        print '\nAdding space charge nodes ...'
        # Make a SC solver
        sizeX = 32
        sizeY = 4
        sizeZ = 4  # Number of longitudinal slices in the 2.5D solver
        calcsbs = SpaceChargeCalcSliceBySlice2D(sizeX,sizeY,sizeZ)
        sc_path_length_min = 0.00000001
        # Add the space charge solver to the lattice as child nodes
        sc_nodes = scLatticeModifications.setSC2p5DAccNodes(Lattice, sc_path_length_min, calcsbs)
        print '  Installed', len(sc_nodes), 'space charge nodes ...'

#-----------------------------------------------------
# Add tune analysis child node
#-----------------------------------------------------
parentnode_number = 97
parentnode = Lattice.getNodes()[parentnode_number]
Twiss_at_parentnode_entrance = Lattice.getNodes()[parentnode_number-1].getParamsDict()
tunes = TeapotTuneAnalysisNode("tune_analysis")

if slicebyslice:
        tunes.assignTwiss(Twiss_at_parentnode_entrance['betax'], Twiss_at_parentnode_entrance['alphax'], Twiss_at_parentnode_entrance['etax'], Twiss_at_parentnode_entrance['etapx'], Twiss_at_parentnode_entrance['betay'], Twiss_at_parentnode_entrance['alphay'])
        addTeapotDiagnosticsNodeAsChild(Lattice, parentnode, tunes)
if frozen:
        tunes.assignTwiss(*[Twiss_at_parentnode_entrance[k] for k in ['betax','alphax','etax','etapx','betay','alphay','etay','etapy']])
        tunes.assignClosedOrbit(*[Twiss_at_parentnode_entrance[k] for k in ['orbitx','orbitpx','orbity','orbitpy']])
        addTeapotDiagnosticsNodeAsChild(Lattice, parentnode, tunes)

#----------------------------------------------------
# Prepare a bunch object to store particle coordinates
#----------------------------------------------------
bunch_tmp = Bunch()
bunch.copyEmptyBunchTo(bunch_tmp)
bunch_tmp.addPartAttr('ParticlePhaseAttributes')

#----------------------------------------------------
# Define twiss analysis and output dictionary
#----------------------------------------------------
bunchtwissanalysis = BunchTwissAnalysis() #Prepare the analysis class that will look at emittances, etc.
get_dpp = lambda b, bta: np.sqrt(bta.getCorrelation(5,5)) / (b.getSyncParticle().gamma()*b.mass()*b.getSyncParticle().beta()**2)
get_bunch_length = lambda b, bta: 4 * np.sqrt(bta.getCorrelation(4,4)) / (speed_of_light*b.getSyncParticle().beta())
get_eps_z = lambda b, bta: 1e9 * 4 * pi * bta.getEmittance(2) / (speed_of_light*b.getSyncParticle().beta())

output = Output_dictionary()
output.addParameter('turn', lambda: turn)
output.addParameter('intensity', lambda: bunchtwissanalysis.getGlobalMacrosize())
output.addParameter('n_mp', lambda: bunchtwissanalysis.getGlobalCount())
output.addParameter('gamma', lambda: bunch.getSyncParticle().gamma())
output.addParameter('mean_x', lambda: bunchtwissanalysis.getAverage(0))
output.addParameter('mean_xp', lambda: bunchtwissanalysis.getAverage(1))
output.addParameter('mean_y', lambda: bunchtwissanalysis.getAverage(2))
output.addParameter('mean_yp', lambda: bunchtwissanalysis.getAverage(3))
output.addParameter('mean_z', lambda: bunchtwissanalysis.getAverage(4))
output.addParameter('mean_dE', lambda: bunchtwissanalysis.getAverage(5))
output.addParameter('epsn_x', lambda: bunchtwissanalysis.getEmittanceNormalized(0))
output.addParameter('epsn_y', lambda: bunchtwissanalysis.getEmittanceNormalized(1))
output.addParameter('eps_z', lambda: get_eps_z(bunch, bunchtwissanalysis))
output.addParameter('bunchlength', lambda: get_bunch_length(bunch, bunchtwissanalysis))
output.addParameter('dpp_rms', lambda: get_dpp(bunch, bunchtwissanalysis))

if frozen or slicebyslice:
        output.addParameter('BE_intensity1', lambda: sc_params1['intensity'])
        output.addParameter('BE_epsn_x1', lambda: sc_params1['epsn_x'])
        output.addParameter('BE_epsn_y1', lambda: sc_params1['epsn_y'])
        output.addParameter('BE_dpp_rms1', lambda: sc_params1['dpp_rms'])

# Define particle output dictionary
#-----------------------------------------------------------------------
particle_output = Particle_output_dictionary()

# Automatically adds particle 0, lets add the rest
# ~ for i in range(1, p['n_macroparticles']):
	# ~ particle_output.AddNewParticle(i)
	
# ~ particle_output.AddNewParticle(1)
# Update for turn -1 (pre tracking)
particle_output.update(bunch, -1)

#----------------------------------------------------
# Do some turns and dump particle information
#----------------------------------------------------
print '\nnow start tracking...'

for turn in range(p['turns_max']):
	Lattice.trackBunch(bunch, paramsDict)
        LinearRestoringForce(bunch, s['RestoringForce'])
        
	bunchtwissanalysis.analyzeBunch(bunch)  # analyze twiss and emittance	

	if turn in p['turns_print']:
		saveBunchAsMatfile(bunch, "output/mainbunch_%s"%(str(turn).zfill(6)))
		saveBunchAsMatfile(lostbunch, "lost/lostbunch_%s"%(str(turn).zfill(6)))

	output.save_to_matfile('output')
	output.update()
	particle_output.update(bunch, turn)
	
	# For last turn output particle dictionary and/or make plots
	if turn == (p['turns_max']-1):
		for i in range(0, p['n_macroparticles']):
			particle_output.print_particle(i)
					
		particle_output.print_all_particles()
		
