#!/bin/csh


nmrPipe -in test.fid \
| nmrPipe -fn EM -lb 5.0    \
| nmrPipe -fn ZF -auto                      \
| nmrPipe -fn FT                            \
| nmrPipe -fn PS -p0 -45.668 -p1 0.0 -di      \
   -out test.ft2 -verb -ov
