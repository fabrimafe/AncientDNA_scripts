#!/bin/bash
source ~/Dropbox/Vindija/scripts/myfunctions.sh --source-only
cmd=$@
echo $cmd
eval $cmd