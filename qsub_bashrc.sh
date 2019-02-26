#!/bin/bash
source ~/.bashrc
source ~/Dropbox/Vindija/scripts/myfunctions.sh --source-only
cmd=$@
echo $cmd
eval $cmd
