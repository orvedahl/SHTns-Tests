#!/bin/sh -p

Ncores=28

output=Multiple-Runs-Timing/Rayleigh

l=255; Nr=128;
for i in $(seq 1 $Ncores); do
  ./xmain.Rayleigh --lmax ${l} --nfields 20 --nr ${Nr} --nloops 10 \
            --output ${output}/timing_grid_${i}.out --n-threads 1 --min-walltime 35.0 \
            --max-walltime 3600.0 > ${output}/out_${i}.out 2>&1 &
done

sleep 35s
