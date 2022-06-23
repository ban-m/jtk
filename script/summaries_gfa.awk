#!/bin/bash awk
($1 ~ /^S/ && $0 ~ /cp:i:1/) { num1 +=1 ; len1 += $3 }
($1 ~ /^S/ && $0 ~ /cp:i:2/) { num2 +=1 ; len2 += $3}
END{
    print "1", num1,len1 / 1000
    print "2", num2,len2 / 1000
    print "total", num1+ num2, (len1 + len2 ) / 1000
}