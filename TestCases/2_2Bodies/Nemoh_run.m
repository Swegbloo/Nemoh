function Nemoh_run()


w = 0.01:0.1:4;
dir = 0;
depth = 0;

[A,B,Fe] = Nemoh('C:\Users\Swegbloo\Documents\GitHub\Nemoh\TestCases\2_2Bodies',1,0);
save("A,B,Fe");
