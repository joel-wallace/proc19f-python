nmrPipe -in testH.fid \
| nmrPipe -fn EM -lb 0.3 -c 0.5  \   # Less line broadening (~0.3 Hz typical for 1H)
| nmrPipe -fn ZF -auto            \   # Zero-filling (same)
| nmrPipe -fn FT -auto            \   # Fourier transform (same)
| nmrPipe -fn PS -p0 0 -p1 0.00 -di -verb  \  # Adjust phase manually
| nmrPipe -fn BASE -nw 10 -nl 0ppm 12ppm  \  # Baseline correction for 1H range
| nmrPipe -fn EXT -x1 0ppm -xn 15ppm -sw  \  # Extract 1H region (0-12 ppm)
   -ov -out test_1H.ft1

