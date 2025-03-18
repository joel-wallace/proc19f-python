#!/bin/csh

bruk2pipe -verb -in ./fid \
  -bad 0.0 -ext -aswap -AMX -decim 2000 -dspfvs 21 -grpdly 76  \
  -xN             32768  \
  -xT             16384  \
  -xMODE            DQD  \
  -xSW        10000.000  \
  -xOBS         500.182  \
  -xCAR           4.771  \
  -xLAB              1H  \
  -ndim               1  \
| nmrPipe -fn MULT -c 4.88281e-01 \
  -out ./test.fid -ov

sleep 5
