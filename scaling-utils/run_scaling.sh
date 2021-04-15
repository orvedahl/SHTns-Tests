#!/bin/sh -p

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# how to use this script:
#     1) generate a test directory
#     2) build the code and setup the input file (with no output)
#     3) setup the PBS/Slurm jobscript to request enough resources
#     4) have the jobscript call this script as:
#            ...
#            bash run_scaling.sh <executable> <mpiexec_command> <output_basename>
#            ...
#     For example, run on Pleiades using their MPI module:
#           bash run_scaling.sh ./xrayleigh mpiexec_mpt scaling_rayleigh
#
#     Run on Pleiades using their MPI module, but with different executable:
#           bash run_scaling.sh ./xrayleigh.SHTns mpiexec_mpt scaling_SHTns
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_it() {
   _exe=$1
   _lmax=$2
   _nr=$3
   _niter=$4
   _out=$5
   _nc=$6
   _mpi_cmd=$7
   ${_mpi_cmd} -np ${_nc} ${_exe} -nr ${_nr} -lmax ${_lmax} -niter ${_niter} > ${_out} 2>&1
}

exe=$1
mpi_cmd=$2
out=$3

# define values for each scaling run: l_max, N_r, N_cores, N_iterations
lvals=( 127 255 511 1023 )
Nrs=( 64 128 )
Ncores=( 8 16 32 64 128 256 512 1024 )
Niter=50

# run each case
for l in "${lvals[@]}"; do
  for Nr in "${Nrs[@]}"; do
    for Nc in "${Ncores[@]}"; do
      run_it ${exe} ${l} ${Nr} ${Niter} ${out}_lmax${l}_Nr${Nr}_Nc${Nc}.out ${Nc} ${mpi_cmd}
      sleep 10s
    done
  done
done

