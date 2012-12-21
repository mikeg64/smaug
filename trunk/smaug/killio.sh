#!/bin/bash

if [ $# -gt 0 ]; then
  IOME_SIMNAME=$1
else
  IOME_SIMNAME="mysim"
fi


INPUT=`cat ${IOME_SIMNAME}0_port.txt`
IOME_WSPORT=$(echo $INPUT | cut -d' ' -f1 )
echo killing iome on port  $IOME_WSPORT



iogs exitiome  0 $IOME_WSPORT
