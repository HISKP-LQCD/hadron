#!/bin/bash
# $Id$

ls outprc* | sort -r -k 3 -t "." | xargs grep "^ [ 0-2][ 0-9]  ." | grep -v "^ 21" | sed -e 's/://g' | sed -e 's/outpr.*\...\.//g' | awk '{print $2 " " $3 " " $4 " " $5 " " $6 " " $1'} > pion.dat
