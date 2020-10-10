# Baseline algorithm of SPO2_HR

# formula:
R = Red(ac)/Red(dc) / Ir(ac)/Ir(dc).  ac: AC component; dc: DC component;

Spo2 = A*R*R + B * R + C (polynomial)

By detecting  peaks of PPG cycle and corresponding AC/DC of red/infra-red signal, the ratio for the SPO2 is computed.
As a reference, this baseline is original code. If you want to use it in your project, you should polish it again.

step 1： update the Freq: FS in algorithm.h

step 2: update the your own R curve in the algorithm.h

step 3： put the files into the visual studio
