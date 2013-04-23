x=(-100:1:100);
y=[];
for i=1:1:length(x)
    y(i)=foo(x(i));
end
plot(x,y);
x0=[50];
options = optimset('Algorithm','sqp');
[x,fval] = fmincon(@foo,x0,[],[],[[10]],[[100]],[],[],[],options);
x 
fval