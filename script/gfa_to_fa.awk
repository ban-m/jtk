#!/usr/bin/awk -f
($1 ~ /^S/){print ">" $2 "\n" $4}