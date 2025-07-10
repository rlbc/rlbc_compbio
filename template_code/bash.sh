#!/bin/bash
# Program.sh
# 
# Program description goes here
#
# Version 0.1: First version
#

#FLAGS

USAGE_MSG="
Usage: $(basename "$0") [OPTIONS]

OPTIONS:
- -t THREADS         Numbers of threads to use (defaults to 1)
- -V                 Version
- -h                 This message
"

# Exits if no options
if [ $# -eq 0 ]; then
  echo "$USAGE_MSG"
  exit 1
fi

# Treatment of command line options
while getopts ":hVt:" options
do
  case "$options" in
    t) THREADS="$OPTARG" ;;
  
    h)
      echo "$USAGE_MSG"
      exit 0
    ;;
  
    V)
      echo -n $(basename "$0")
      # Extract version from program heading
      grep '^# Version ' "$0" | tail -1 | cut -d : -f 1 | tr -d \#
      exit 0
    ;;

    \?)
      echo Invalid option: $OPTARG
      exit 1
    ;;

    :)
      echo Missing argument for: $OPTARG
      exit 1
    ;;
  esac
done

MYTIME=$(date +%Y%m%d%-H%M%S)

# From here push logs to specific file
exec 1>>$(basename "$0")_${MYTIME}.log
exec 2>>$(basename "$0")_${MYTIME}_2.log

##############
# CODE.......#
##############
