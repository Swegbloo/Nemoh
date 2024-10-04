function Nemoh_run()
clear all;
[t,w,A,B,Fe] = Nemoh('C:\Users\Swegbloo\Documents\GitHub\Nemoh\sdc_cy\1_Cylinder',1,0);
save('Nemohresult_test','A','B','Fe','w','t');
