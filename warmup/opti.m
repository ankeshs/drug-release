x0 = [-1,1];     % Make a starting guess at the solution
options = optimset('Algorithm','active-set');
[x,fval] = fmincon(@objfun1,x0,[],[],[],[],[],[],@confun1,options);
x
fval