#!/bin/bash
#-------------------------------------------------------------------------------
# WABBIT unit test
# This file contains one specific unit test, and it is called by unittest.sh
#-------------------------------------------------------------------------------
# what parameter file
dir="./TESTING/acm/acm_cyl/"
params=${dir}"acm_test.ini"
happy=0
sad=0
mpi_command="mpiexec -n 4"
echo "testing artificial compressibility"

# list of prefixes the test generates
prefixes=(Ux Uy p mask vor)
# list of possible times (no need to actually have them)
times=(000000002000)

# run actual test
${mpi_command} ./wabbit 2D ${params} --memory=2GB

echo "============================"
echo "run done, analyzing data now"
echo "============================"

# loop over all HDF5 files and generate keyvalues using wabbit
for p in ${prefixes[@]}
do
  for t in ${times[@]}
  do
    echo "--------------------------------------------------------------------"
    # *.h5 file coming out of the code
    file=${p}"_"${t}".h5"
    # will be transformed into this *.key file
    keyfile=${p}"_"${t}".key"
    # which we will compare to this *.ref file
    reffile=${dir}${p}"_"${t}".ref"

    if [ -f $file ]; then
        # get four characteristic values describing the field
        ${mpi_command} ./wabbit-post 2D --keyvalues ${file}
        # and compare them to the ones stored
        if [ -f $reffile ]; then
            ${mpi_command} ./wabbit-post 2D --compare-keys $keyfile $reffile
            result=$(cat return); rm return
            if [ $result == "0" ]; then
              echo -e "\033[92m :) Happy, this looks okay!\033[0m" $keyfile $reffile
              happy=$((happy+1))
            else
              echo -e "\033[31m:[ Sad, this is failed!\033[0m" $keyfile $reffile
              sad=$((sad+1))
            fi
        else
            sad=$((sad+1))
            echo -e "\033[31m:[ Sad: Reference file not found\033[0m"
        fi
    else
        sad=$((sad+1))
        echo -e "\033[31m:[ Sad: output file not found\033[0m"
    fi
    echo " "
    echo " "

  done
done


echo -e "\t \033[92m happy tests:\033[0m \t" $happy
echo -e "\t \033[31m sad tests:\033[0m \t" $sad

#-------------------------------------------------------------------------------
#                               RETURN
#-------------------------------------------------------------------------------
if [ $sad == 0 ]
then
  exit 0
else
  exit 999
fi