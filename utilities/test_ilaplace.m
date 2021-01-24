close all;
clear
syms s
t1 = 0:0.01:10;
x = 5*exp(-1);
y1 = 5*(exp(-2*t1)-exp(-t1));
y2sym = ilaplace(5*s/(s+2)/(s+1)-1/(s+2)*5);
y2 = vpa(subs(y2sym,t1),3);
figure; plot(t1,y1,t1,y2);legend('y1','y2');

ilaplace(5/(s+2))

