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
   echo "  --grid         Run grid of tests"
   echo "  --high-res     Run test with larger lmax & Nr"
   echo "  --more-loops   Increase minimum number of loops to use"
   echo ""
}

parse_cmd() {
   exe=0
   output=0
   high_res=0
   Nloops=10
   run_grid=0

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
           --grid )
              run_grid=1
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
      --output "${_out}" --n-threads ${_nthread} --min-walltime 35.0 --max-walltime 1800.0
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
   Nrs=( 4 8 16 32 64 )
fi

if [ "$run_grid" = "0" ]; then
   # run timing tests
   for l in "${lvals[@]}"; do
      run_it ${exe} ${l} 16 10 1 ${output}/timing_fixed_Nrhs.out ${Nloops}
   done

   for Nr in "${Nrs[@]}"; do
      run_it ${exe} 63 ${Nr} 10 1 ${output}/timing_fixed_lmax.out ${Nloops}
   done
else
   # run grid
   for l in "${lvals[@]}"; do
      for Nr in "${Nrs[@]}"; do
         run_it ${exe} ${l} ${Nr} 20 1 ${output}/timing_grid.out ${Nloops}
      done
   done
fi

