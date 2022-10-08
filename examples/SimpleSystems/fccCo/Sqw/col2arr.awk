#!/bin/bash
awk 'BEGIN {RS="";FS="\n";ORS=" "}{ dim=3; for (i=0; i<= NF; i++) {
                           print $(i+1)
                           if ((i%diml )==diml-1) printf "\n"
                                                }}' diml=$2  $1
