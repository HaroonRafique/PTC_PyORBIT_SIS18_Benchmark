#!/bin/bash
# submit job using: condor_submit htcond.sub

# some MPI related stuff
export MPI=/afs/cern.ch/project/lhcnap/public/public_html/Source/space_charge/Codes/mpich2-1.4.1p1
export MPIBIN=$MPI/bin

# get condor variables
OrigIwd=`grep "^OrigIwd" .job.ad | cut -d'"' -f2`
ClusterId=`grep "^ClusterId" .job.ad | tr -dc '0-9'`
ProcId=`grep "^ProcId" .job.ad | tr -dc '0-9'`

# retreiver_file=${OrigIwd}"_retreive_summary.sh"
# echo "#!/bin/bash" > $retreiver_file
# echo "condor_ssh_to_job -ssh scp ${ClusterId}.${ProcId} remote:output/output.mat ${OrigIwd}" >> ${retreiver_file}
# chmod u+x $retreiver_file

# cd to original working directory
cd ${OrigIwd}

input_dir='input'
output_dir='output'
mkdir -p $input_dir
mkdir -p $output_dir
simulation_info_file="${output_dir}/simulation_info_${ClusterId}.${ProcId}.txt"
echo "PyOrbit path:  `readlink -f ${ORBIT_ROOT}`" > ${simulation_info_file}
echo "using $1 CPUs" >> ${simulation_info_file}
echo "****************************************" >> ${simulation_info_file}

# copy all files from initial directory
# cp -r ${OrigIwd}* .

# run the job
tstart=$(date +%s)
nnodes=$1
$MPIBIN/mpirun -np $nnodes ${ORBIT_ROOT}/bin/pyORBIT pyOrbit.py
tend=$(date +%s)
dt=$(($tend - $tstart))
echo 'total simulation time (s): ' $dt >> ${simulation_info_file}

# copy files back to initial directory
# cp -r ${input_dir} ${output_dir} ptc ${OrigIwd}

rm Maxwellian_bend_for_ptc.txt
rm junk.txt
