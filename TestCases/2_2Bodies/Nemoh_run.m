function Nemoh_run()


w = 0.01:0.1:4;
dir = 0;
depth = 0;

[~,~,A,B,Fe] = Nemoh(w, dir, depth);
