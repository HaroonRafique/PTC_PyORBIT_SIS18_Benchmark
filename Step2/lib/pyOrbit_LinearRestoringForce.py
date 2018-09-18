#-----------------------------------------------------------------------
# Function for linear restoring force, replacing normal RF cavity action
# 18.09.2018: Created by Haroon Rafique, CERN BE-ABP-HSI 
#-----------------------------------------------------------------------

import math
import orbit_mpi

from orbit_mpi import mpi_datatype, mpi_op
from bunch import Bunch

def LinearRestoringForce(bunch, force):

	rank = 0
	numprocs = 1

	mpi_init = orbit_mpi.MPI_Initialized()
	comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD

	if(mpi_init):
		rank = orbit_mpi.MPI_Comm_rank(comm)
		numprocs = orbit_mpi.MPI_Comm_size(comm)

	nparts_arr_local = []
	for i in range(numprocs):
		nparts_arr_local.append(0)
			
	nparts_arr_local[rank] = bunch.getSize()
	data_type = mpi_datatype.MPI_INT
	op = mpi_op.MPI_SUM

	nparts_arr = orbit_mpi.MPI_Allreduce(nparts_arr_local,data_type,op,comm)

	for i in range(bunch.getSize()):
		en = bunch.dE(i)

		en = en + bunch.z(i) * force
		
		bunch.dE(i,en)
		
	return
