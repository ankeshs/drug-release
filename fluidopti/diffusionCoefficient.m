function [ D Ap Bp Cp Dp ] = diffusionCoefficient( D01,Q1,X,V1,V2,K21mTg2,K22mTg2,E,K11dY, K12dY, T,w1,w2 )
%DIFFUSIONCOEFFICIENT Summary of this function goes here
%   Detailed explanation goes here
    t1 = (1-Q1^2)*(1-2*X*Q1);
    Ap=E*V2;
    Bp=V1-E*V2;
    Cp=K12dY*(K22mTg2+T);
    Dp=K11dY*(K21mTg2+T)-K12dY*(K22mTg2+T);
    t2 = -(Ap+Bp*w1);
    t3 = Cp+Dp*w1;
    D = D01*t1*exp(t2/t3);
end

