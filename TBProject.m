function fval = TBProject(t,y)
% Function for TBmodel system
% This set of ODEs 
S = y(1);
L1 = y(2);
I = y(3);
L2 = y(4);
R = y(5);

% Define Constants
b  = 100;
miu = 1/52;
d = 12;
phi = 0.05;
w = 0.0002;
wr = 0.00002;
rho = 0.25;
rhor = 0.25;
tao0 = 2;
tao1 = 2;
tao2 = 1;
N = 30000;
e1 = 0.5;
e2 = 0.5;
% Will decide later
U1 = 0.1;
U2 = 0.1;

% T = 5;
% W1 = 500;
% W2 = 50;


% Define dy/dtd

fval(1,1) = miu*N - b/N*I*S - miu*S;
fval(2,1) = b/N*I*(S + rho*L2 + rhor*R) - (d + tao1 + miu)* L1;
fval(3,1) = phi*d*L1 + w*L2 + wr*R - (tao0 + e1*U1 + miu)*I;
fval(4,1) = (1 - phi)*d*L1 - rho*b/N*I*L2 - (w + e2*U2 + tao2 + miu)*L2;
fval(5,1) = (tao0 + e1*U1)*I + tao1*L1 + (tao2 + e2*U2)*L2 - rhor*b/N*I*R-(wr + miu)*R;







