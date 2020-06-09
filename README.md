# HLA class

Author: Bansho Masutani

Mail: ban-m@g.ecc.u-tokyo.ac.jp

# Desc

## Implementation details


### Unit selection

1. Sort reads in descending order w.r.t the length of the read.
2. Until we collect sufficient number of chunks or there is no read, do:
   1. Pop the longest read.
   2. Take the L-length chunk every K-base skipping interval



### Unit encoding 

# Memos

/data/hacone/randseq/