# $Id$

ls onlinemeas* | sort -k 2 -t "." | xargs grep "^[0-6][ 0-9] ." | grep -v "^21" | sed -e 's/:/ /g' | sed -e 's/onlinemeas\.//g' | awk '{print $2 " " $3 " " $4 " " $5 " " $6 " " $1'} > piononline.dat
