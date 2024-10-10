clc
clear all

nfreq = 50;
w = linspace(2.73, 8.86, nfreq);
dir = 0;
depth = 1;
[A,B,Fe] = Nemoh(w, dir, depth);
save('Nemohresult_test','A','B','Fe','w');