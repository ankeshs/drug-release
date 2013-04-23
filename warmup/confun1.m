function [c, ceq] = confun1(x)
% Nonlinear inequality constraints
c = [1.5 + x(1)*x(2) - x(1) - x(2);     
     -x(1)*x(2) - 10];
% Nonlinear equality constraints
ceq = [];