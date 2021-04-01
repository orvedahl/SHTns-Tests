#!/bin/sh -p

# function to run both SHTns and Rayleigh
run_it() {
   lmax=$1
   nr=$2
   nf=$3
   nloops=$4
   nthread=$5
   out=$6
   ./xmain.Linux.intel.Rayleigh.exe --run-timing 1 --nr ${nr} --lmax ${lmax} --nfields ${nf} --nloops ${nloops} --output "${out}" --n-threads ${nthread} --walltime 600.0
   #./xmain.Linux.intel.SHTns.exe --run-timing 1 --nr ${nr} --lmax ${lmax} --nfields ${nf} --nloops ${nloops} --output "${out}" --n-threads ${nthread} --walltime 600.0
}


lvals=( 15 31 63 127 255 511 1023)
for l in "${lvals[@]}"; do
   run_it ${l} 16 10 30 1 timing_fixed_Nrhs_320.Rayleigh.out
done

Nrs=( 2 4 8 16 32 64 128 256)
for Nr in "${Nrs[@]}"; do
   run_it 63 ${Nr} 10 30 1 timing_fixed_lmax_63.Rayleigh.out
done

Nfs=( 2 4 8 16 32 64 128 256)
for Nf in "${Nfs[@]}"; do
   run_it 63 16 ${Nf} 30 1 timing_fixed_lmax63_Nr32.Rayleigh.out
done


###############################################################################
# loop over a bunch of variables to get timing
#lvals=( 15 31 63 127 )
#Nrs=( 16 32 64 128 )
#Nfs=( 5 10 20 30 )
#Nloops=( 3 10 30 100 300 )
#Nthreads=( 1 2 4 8 )
#
#for l in "${lvals[@]}"; do
#  for Nr in "${Nrs[@]}"; do
#     for f in "${Nfs[@]}"; do
#        for lo in "${Nloops[@]}"; do
#           for n in "${Nthreads[@]}"; do
#              run_it ${l} ${Nr} ${f} ${lo} ${n} timing_nthreads_${n}.out
#           done
#        done
#     done
#  done
#done
#
## now do the Nloop = 1000
#lo=1000
#for l in "${lvals[@]}"; do
#  for Nr in "${Nrs[@]}"; do
#     for f in "${Nfs[@]}"; do
#        for n in "${Nthreads[@]}"; do
#           run_it ${l} ${Nr} ${f} ${lo} ${n} timing_nthreads_${n}.out
#        done
#     done
#  done
#done
