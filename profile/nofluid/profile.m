Di=1e-14;
Do=1e-15;
Kr=0.4;
c0=100;
ar=0.4;
cd=0.2;
tf=15000;
dt=5;
st=50;
M=condep(c0, 0, 0, 1e-11, Kr, tf, dt);
n=M.result.numerical('int1').getReal;
t1=0:20:1000;
t2=1000:1000:15000;
t=[t1 t2];
t
n
plot(t,n);