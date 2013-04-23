function [c ceq] = constraint(X)

D01 = X(1);		% [1e-8 1e+8]
Q1 = 0;			
V1 = 0.75;
V2 = X(2);		% [0.6 1.3]
K21mTg2 = X(3);		% [-400 -50]
K22mTg2 = X(4);		% [-100 -10]
E = X(5);		% [0.1 10]
K = 2.59;		
K11dY = 2.72e-5;
K12dY = X(6)*1e-4;	% [1 10]
Kr = X(7);		% [0.05 0.4]
hp = X(8);		% [1e-3 0.35]
Co = 100;
rho = X(9);		% [900 1200]
T = 310;

[ D1 Ap Bp Cp Dp ] = diffusionCoefficient( D01, Q1, K, V1, V2, K21mTg2, K22mTg2, E, K11dY, K12dY, T, 0.01, 0.99 );
[ D2 Ap Bp Cp Dp ] = diffusionCoefficient( D01, Q1, K, V1, V2, K21mTg2, K22mTg2, E, K11dY, K12dY, T, 0.1, 0.9 );
[ D3 Ap Bp Cp Dp ] = diffusionCoefficient( D01, Q1, K, V1, V2, K21mTg2, K22mTg2, E, K11dY, K12dY, T, 0.2, 0.8 );

% Constraints
% 1e-15 <= D1, D2, D3 <= 1e-8
% 1e-15 - Di <= 0 and Di - 1e-8 <= 0 

c = [1e-15-D1 ; D1-1e-8 ; 1e-15-D2 ; D2-1e-8 ; 1e-15-D3 ; D3-1e-8 ; ];
ceq = [];

end