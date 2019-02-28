#!/bin/bash
awk -v OFS=' ' -v myrescale=$THETAFACTORNEW '{ if ($1=="T"){ print $1,$2*myrescale} else print}'
