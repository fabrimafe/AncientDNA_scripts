#!/bin/bash
#argument is $THETAFACTORNEW
awk -v OFS=' ' -v myrescale=$1 '{ if ($1=="T"){ print $1,$2*myrescale} else print}'
