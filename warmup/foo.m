function fval = foo( x )
xa=-1000:1:1000;
ya=[];
for i=1:1:length(xa)
    ya(i)=5*xa(i)^2-0.005*xa(i)^3+100*xa(i);
end
d=2000;
di=0;
for i=1:1:length(xa)
    t=abs(x-xa(i));
    if(t<d)
        d=t;
        di=i;
    end
end
si=di;
ei=di+1;
if(x < xa(di)) 
    si=di-1;
    ei=di;
end
fval=ya(si)+(ya(ei)-ya(si))*(x-xa(si))/(xa(ei)-xa(si));
end

