Di=1e-14;
Do=1e-15;
Kr=0.4;
c0=100;
ar=0.4;
cd=0.2;
tf=6000;
dt=20;
st=50;
Df=3e-8;
Kp=0.1;


t=[0:st:tf];
d=[1e-16 5e-16 1e-15 5e-15 1e-14];
clear ln;
ln=[];
for i=1:length(d)
M=polyff(Di, d(i), Df, Kr, Kp, c0, ar, cd, tf, dt, st);
n=M.result.numerical('int1').getReal;
p=plot(t,n);
hold all;
set(p,'LineWidth',2);
ln=strvcat(ln,strcat('Do = ', num2str(d(i))));
end
title('Release Profile');
xlabel('Time (s)');
ylabel('Drug release (mol/s)');
legend(cellstr(ln));