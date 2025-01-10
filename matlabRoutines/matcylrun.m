clc
clear all

addpath('C:\Users\swago\OneDrive\Documents\GitHub\Nemoh\matlabRoutines\cylinder');

nfreq = 2;
w = linspace(2.73, 8.86, nfreq);
dir = 0;
depth = 0;
[A,B,Fe] = Nemoh(w, dir, depth);
save('Nemohresult_test','A','B','Fe','w');