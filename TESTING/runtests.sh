#!/bin/bash

#--------------------------------
# WABBIT run unit tests
#--------------------------------
rootdir="../"
cd $rootdir
echo "   "
echo -e "\t \033[4m WABBIT: run all existing unit tests \033[0m" 
echo "   "

if [ -z "$nprocs" ]; then
    nprocs=4
fi

if [ -z "$mpi_command" ]; then
    export mpi_command="nice mpiexec -n ${nprocs}"
fi

fail_color=$'\033[31;1m'
pass_color=$'\033[92;1m'
end_color=$'\033[0m'

tests=("---navier-stokes---" "TESTING/navier_stokes/pressure_blob/pressure_blob.sh" "---acm---" "TESTING/acm/acm_cyl/acm_cylinder.sh")
happy_sum=0
sad_sum=0
numtests=0

echo "employed command for parallel exec: " $mpi_command
echo "to modify the command, set \$nprocs, \$mpi_command in shell"
echo "   "

T="$(date +%s)"

for ts in ${tests[@]}
do
    if [[ $ts == "---"* ]]; then
    	echo $ts
    else
	logfile=${ts%%.sh}.log
        rm -f $logfile
	touch $logfile

	./${ts} > $logfile
	if [ $? == 0 ]; then
	    happy_sum=$((happy_sum+1))
	    summary[$numtests]=0
	else
	    sad_sum=$((sad_sum+1))
	    summary[$numtests]=1
	fi
	numtests=$((numtests+1))
	rm -f *.key *.h5 *.t *.dat
    fi
done

echo
T="$(($(date +%s)-T))"
echo "Time used in tests: ${T} seconds"

echo " "
echo "All in all we have: "
echo " "

numtests=0
for ts in ${tests[@]}
do
    if [[ $ts != "---"* ]]; then
	if [ ${summary[$numtests]} == 0 ]; then
	    printf "%-80s %s \n" ${ts} "$pass_color ok $end_color"
	else
            printf "%-80s %s \n" ${ts} "$fail_color X $end_color"
	fi
	numtests=$((numtests+1))
    fi
done
echo " "

echo -e "\t $pass_color sum happy tests: $end_color \t" $happy_sum
echo -e "\t $fail_color sum sad tests: $end_color \t" $sad_sum