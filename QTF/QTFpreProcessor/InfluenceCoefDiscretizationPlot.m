%plotting discretisation for second term of influence coefficient
%In NEMOH1
close all;clear all;
NI=328;
NJ=46;
Zmin=-16;
Z=zeros(NJ,1);
X=zeros(NI,1);
for J=1:NJ
 Z(J)=-min(min(10^(J/5-6),10^(J/8-4.5)),-Zmin);
end
for I=2:NI
 % if I<32
  X(I)=min(10^((I-1)/5-6),4/3+abs((I-32))/3);
 % else
%  X(I)=4/3+abs((I-32))/3;    
 % end
end
figure;
subplot(2,1,1)
plot(X)
ylabel('X')
xlabel('I')
subplot(2,1,2)
plot(Z)
ylabel('Z')
xlabel('J')
%In NEMOH2
NI=700;
NJ=130;
Zmin=-251;
Xmax=100;
Z=zeros(NJ,1);
X=zeros(NI,1);
for J=1:NJ
 Z(J)=-min(10^(J/10-10),-Zmin);
end
for I=2:NI
%  if I<32
  X(I)=min(min(10^((I-1)/10-8),abs((I-75))/6),Xmax);
%  else
%  X(I)=4/3+abs((I-32))/3;    
%  end
end

figure;
subplot(2,1,1)
plot(X)
ylabel('X')
xlabel('I')
subplot(2,1,2)
plot(Z)
ylabel('Z')
xlabel('J')