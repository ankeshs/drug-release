function t = objective(X)
%WINOBJ Summary of this function goes here
%   Detailed explanation goes here
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

[ D Ap Bp Cp Dp ] = diffusionCoefficient( D01, Q1, K, V1, V2, K21mTg2, K22mTg2, E, K11dY, K12dY, T, 0.05, 0.95 );

Rmec=1e-18;
Rtox=1e-17;

t=-1*thWindow([Co Ap Bp Cp Dp D01 Kr rho hp 15000 10], ['conc' 'A' 'B' 'C' 'D' 'D01' 'poly partition' 'density' 'thickness' 'duration' 'time step'], Rmec, Rtox);

% Overall variable limits for optimization:
  % [1e-8 0.6 -400 -100 0.1 1 0.05 1e-3 900] to [1e+8 1.3 -50 -10 10 10 0.4 0.35 1200]

% Generate DB
fh = fopen ('runData2.dat','a');
fprintf(fh, '%f %f %f %f %f %f %f %f %f %f\n', X(1), X(2), X(3), X(4), X(5), X(6), X(7), X(8), X(9), -1*t);
fclose(fh);

end

