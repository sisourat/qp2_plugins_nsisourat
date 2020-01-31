#!/bin/bash

while read x; do             # read current line from file, see last line
  if [[ $x == "" ]]; then    # empty line?
    cat f7spd                # print KBJ function
    echo ""
  else
    echo "$x"	  
  fi
done
