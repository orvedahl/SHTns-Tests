#!/bin/sh -p

usage() {
   echo ""
   echo "Run the timing test for the given executable"
   echo ""
   echo "Usage:"
   echo "    run_timing_test.sh [options] --exe=<exe> --output=<directory>"
   echo ""
   echo "    <exe> : required, the executable to run"
   echo "    <directory> : required, where to store the timing results"
   echo ""
   echo "Options:"
   echo "  --high-res     Run test with larger lmax & Nr"
   echo "  --more-loops   Run tests with slightly more loops"
   echo ""
}

parse_cmd() {
   exe=0
   output=0
   high_res=0
   Nloops=10

   # actual parsing of command line:
   for i in "$@"; do
       case $i in
           -h | --h | --help | -help )
              usage
              exit
              ;;
           --high-res )
              high_res=1
              ;;
           --more-loops )
              Nloops=30
              ;;
           --exe=* )
              exe="${i#*=}"  # parse out the "=" and keep the specified option
              ;;
           --output=* )
              output="${i#*=}"  # parse out the "=" and keep the specified option
              ;;
           * )
              echo
              echo "===================="
              echo "Unrecognized option"
              echo "===================="
              usage
              exit
              ;;
      esac
   done
}

run_it() {
   _exe=$1
   _lmax=$2
   _nr=$3
   _nf=$4
   _nthread=$5
   _out=$6
   _nloops=$7
   ./$1 --run-timing 1 \
      --nr ${_nr} --lmax ${_lmax} --nfields ${_nf} --nloops ${_nloops} \
      --output "${_out}" --n-threads ${_nthread} --walltime 600.0
}


parse_cmd "$@"


if [ "$exe" = "0" ]; then
   echo
   echo "===================="
   echo "Must provide exe"
   echo "===================="
   usage
   exit
fi
if [ "$output" = "0" ]; then
   echo
   echo "===================="
   echo "Must provide output"
   echo "===================="
   usage
   exit
fi

if [ "$high_res" = "1" ]; then
   lvals=( 15 31 63 127 255 511 1023 )
   Nrs=( 2 4 8 16 32 64 128 256 )
else
   lvals=( 15 31 63 127 )
   Nrs=( 2 4 8 16 32 64 )
fi

# run timing tests
for l in "${lvals[@]}"; do
   run_it ${exe} ${l} 16 10 1 ${output}/timing_fixed_Nrhs.out ${Nloops}
done

for Nr in "${Nrs[@]}"; do
   run_it ${exe} 63 ${Nr} 10 1 ${output}/timing_fixed_lmax.out ${Nloops}
done

