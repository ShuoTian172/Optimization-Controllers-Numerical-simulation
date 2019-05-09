toms t
p = tomPhase('p', t, 0, 5, 50);
setPhase(p);
tomStates S2 L12 I2 L22 R2
tomControls U22

% constants
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
N = 30000;

W1 = 500;
W2 = 50;
U12=0;

% Initial guess
% Note: The guess for t_f must appear in the list before expression involving t.
x0 = {collocate(U22==0.001)};

% Box constraints
cbox = {
0 <= collocate(U22) <= 1};

% Boundary constraints
cbnd = {initial({S2 == N*76/120; L12 == N*36/120
I2 == N*5/120; L22 == N*5/120; R2 == N*1/120})};

% ODEs and path constraints
ceq = collocate({
dot(S2) == miu*N - b/N*I2*S2 - miu*S2;
dot(L12) == b/N*I2*(S2 + rho*L22 + rhor*R2) - (d + tao1 + miu)* L12;
dot(I2) == phi*d*L12 + w*L22 + wr*R2 - (tao0 + e1*U12 + miu)*I2;
dot(L22) == (1 - phi)*d*L12 - rho*b/N*I2*L22 - (w + e2*U22 + tao2 + miu)*L22;
dot(R2) == (tao0 + e1*U12)*I2 + tao1*L12 + (tao2 + e2*U22)*L22 - rhor*b/N*I2*R2-(wr + miu)*R2});


% Objective
objective = integrate(I2 + L22 + 0.5*W1*U12^2 + 0.5*W2*U22^2);

% Solve the problem
options = struct;
solution = ezsolve(objective, {cbox, cbnd, ceq}, x0, options);
t = subs(collocate(t),solution);
S2 = subs(collocate(S2),solution);
L12 = subs(collocate(L12),solution);
I2 = subs(collocate(I2),solution);
L22 = subs(collocate(L22),solution);
R2 = subs(collocate(R2),solution);
% U12 = subs(collocate(U12),solution);
U22 = subs(collocate(U22),solution);




subplot(3,2,1)
plot(t,(I+L2)/N,'g--*',t,(I2+L22)/N,'b-');
legend('Two controls','Control U2 only');
title('(A) (I^*+L_2^*)/N');
xlabel('Time t'); 
ylabel('(I+L2)/N');

subplot(3,2,2)
plot(t,S/N,'g--*',t,S2/N,'b-');
legend('Two controls','Control U2 only');
title('(B) S^*/N');
xlabel('Time t'); 
ylabel('S/N');

subplot(3,2,3)
plot(t,L1/N,'g--*',t,L12/N,'b-');
legend('Two controls','Control U2 only');
title('(C) L_1^*/N');
xlabel('Time t'); 
ylabel('L_1/N');

subplot(3,2,4)
plot(t,I/N,'g--*',t,I2/N,'b-');
legend('Two controls','Control U2 only');
title('(D) I^*/N');
xlabel('Time t'); 
ylabel('I/N');

subplot(3,2,5)
plot(t,L2/N,'g--*',t,L22/N,'b-');
legend('Two controls','Control U2 only');
title('(E) L_2^*/N');
xlabel('Time t'); 
ylabel('L_2/N');

subplot(3,2,6)
plot(t,R/N,'g--*',t,R2/N,'b-');
legend('Two controls','Control U2 only');
title('(F) R^*/N');
xlabel('Time t'); 
ylabel('R/N');


