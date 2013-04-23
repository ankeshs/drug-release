Di=1e-14;
Do=1e-15;
Kr=0.4;
c0=100;
ar=0.4;
cd=0.2;
tf=6000;
dt=5;
st=50;

t=[0:st:tf];
d=[1e-16 0.5e-15 1e-15 5e-15 1e-14];

ln=[];
for i=1:length(d)
M=polynf(Di, d(i), Kr, c0, ar, cd, tf, dt, st);
n=M.result.numerical('int1').getReal;
maxn=max(n);
n=n/maxn;
p=plot(t,n);
hold all;
set(p,'LineWidth',2);
ln=[ln; strcat('Di = ',num2str(d(i)))];
end
title('Release Profile');
xlabel('Time (s)');
ylabel('Drug release (mol/s)');
legend(cellstr(ln));