#!/bin/sh -p

# function to run both SHTns and Rayleigh
run_it() {
   lmax=$1
   nr=$2
   nf=$3
   nloops=$4
   nthread=$5
   out=$6
   ./xmain.Linux.gfortran.SHTns.theta.gcc.exe --run-timing 1 --nr ${nr} --lmax ${lmax} --nfields ${nf} --nloops ${nloops} --output "${out}" --n-threads ${nthread} --walltime 600.0
}

output=timing-results/gcc-openmp-enabled-theta-contiguous-Run-with-gcc-onthefly

lvals=( 15 31 63 127 255)
for l in "${lvals[@]}"; do
   run_it ${l} 16 10 30 1 ${output}/timing_fixed_Nrhs_320.out
done

Nrs=( 2 4 8 16 32 64 128)
for Nr in "${Nrs[@]}"; do
   run_it 63 ${Nr} 10 30 1 ${output}/timing_fixed_lmax_63.out
done

#Nfs=( 2 4 8 16 32 64 128 256)
#for Nf in "${Nfs[@]}"; do
#   run_it 63 16 ${Nf} 30 1 ${output}/timing_fixed_lmax63_Nr32.out
#done

