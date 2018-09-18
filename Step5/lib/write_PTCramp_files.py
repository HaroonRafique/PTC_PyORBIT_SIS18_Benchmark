import metaclass
import orbit_mpi

def write_QuadRamp_files(target_file, twissfile, pattern, ptc_source_table):
	comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
	rank = orbit_mpi.MPI_Comm_rank(comm)
	if not rank:
		t = metaclass.twiss(twissfile)
		q_i = [i for i, n in enumerate(t.NAME) if pattern in n]
		with open(target_file, 'w') as fid:
			fid.write('SET ORBIT RAMPING \n')
			for i in q_i:
				fid.write(' ramp\n %s\t"%s"\t%1.9f \n'%(t.NAME[i], ptc_source_table, (t.K1L[i]/t.L[i])/(t.K1L[q_i[0]]/t.L[q_i[0]])))
			fid.write('return')
	orbit_mpi.MPI_Barrier(comm)

def write_SextupoleRamp_files(target_file, pattern, ptc_source_table):
	comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
	rank = orbit_mpi.MPI_Comm_rank(comm)
	if not rank:
		with open(target_file, 'w') as fid:
			fid.write('SET ORBIT RAMPING \n')
			fid.write(' ramp\n %s\t"%s"\t%1.9f \n'%(pattern, ptc_source_table, 1.0))
			fid.write('return')
	orbit_mpi.MPI_Barrier(comm)