#!/bin/bash


help()
{
	echo "This program takes output smaug data and gathers into a single outputfile"
        echo "./smauggather.sh numproc start finish step"

	echo "numproc is the number of processors"
        echo "start is the starting step"
        echo "finish is the final step"
        echo "step is optional and is the interval"

	exit 1
}




#module add libs/cuda/4.0.17
#module add mpi/gcc/openmpi/1.4.4

start=$2
end=$3 
step=1
current=$start



numproc=$1

if test $# -lt 3
then
	help
fi

if test $1 -eq "-help"
then
	help
fi

if test $# -gt 3
then
	step=$4
fi



while [ "$current" -le $end ]
do


        mpirun -np $numproc bin/smaug -o gather $current
        current=`expr $current + $step`

done
