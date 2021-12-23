#!/bin/awk
($1 ~ /^@/){
    print("@" int(NR/4) "," $2 , $4,$5,$6)
    next 
}
{
    print($0)
}
