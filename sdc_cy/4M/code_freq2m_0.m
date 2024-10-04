%clc;
%clear all;
m0 = zeros(1,50);
%f = zeros(1,50);
f = transpose(f);
global v;
fun = @freq2m_0;
 t = 1;
 disp('X');
for i = 1:50
    
    v = f(i)^2;
    x_0 = [0 9]; % interval where the root lie
    
    m0(t) = fzero(fun,x_0);
    t = t+1;
end