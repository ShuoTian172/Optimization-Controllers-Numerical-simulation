toms t
p = tomPhase('p', t, 0, 5, 50);
setPhase(p);
tomStates S0 L10 I0 L20 R0
% tomControls U01 U02

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
U01=0;
U02=0;


% W1 = 500;
% W2 = 50;


% Initial guess
% Note: The guess for t_f must appear in the list before expression involving t.
% x0 = {collocate(U02==0),...
% collocate(U01==0)};
% 
% Box constraints
% cbox = {0 <= collocate(U01) <= 1
% 0 <= collocate(U02) <= 1};

% Boundary constraints
cbnd = {initial({S0 == N*76/120; L10 == N*36/120
I0 == N*5/120; L20 == N*5/120; R0 == N*1/120})};

% ODEs and path constraints
ceq = collocate({
dot(S0) == miu*N - b/N*I0*S0 - miu*S0;
dot(L10) == b/N*I0*(S0 + rho*L20 + rhor*R0) - (d + tao1 + miu)* L10;
dot(I0) == phi*d*L10 + w*L20 + wr*R0 - (tao0 + e1*U01 + miu)*I0;
dot(L20) == (1 - phi)*d*L10 - rho*b/N*I0*L20 - (w + e2*U02 + tao2 + miu)*L20;
dot(R0) == (tao0 + e1*U01)*I0 + tao1*L10 + (tao2 + e2*U02)*L20 - rhor*b/N*I0*R0-(wr + miu)*R0});


% Objective
objective = integrate(I0 + L20+0.5*W1*U01^2 + 0.5*W2*U02^2);

% S0olve the problem
options = struct;
solution = ezsolve(objective, {cbox, cbnd, ceq}, x0, options);
t = subs(collocate(t),solution);
S0 = subs(collocate(S0),solution);
L10 = subs(collocate(L10),solution);
I0 = subs(collocate(I0),solution);
L20 = subs(collocate(L20),solution);
R0 = subs(collocate(R0),solution);
% U01 = subs(collocate(U01),solution);
% U02 = subs(collocate(U02),solution);



subplot(3,2,1)
plot(t,(I+L2)/N,'g--*',t,(I1+L21)/N,'b-',t,(I0+L20)/N,'r.');
legend('Two controls','Control U1 only','Without controls');
title('(A) (I^*+L_2^*)/N');
xlabel('Time t'); 
ylabel('(I+L2)/N');

subplot(3,2,2)
plot(t,S/N,'g--*',t,S1/N,'b-',t,S0/N,'r.');
legend('Two controls','Control U1 only','Without controls');
title('(B) S^*/N');
xlabel('Time t'); 
ylabel('S/N');

subplot(3,2,3)
plot(t,L1/N,'g--*',t,L11/N,'b-',t,L10/N,'r.');
legend('Two controls','Control U1 only','Without controls');
title('(C) L_1^*/N');
xlabel('Time t'); 
ylabel('L_1/N');

subplot(3,2,4)
plot(t,I/N,'g--*',t,I1/N,'b-',t,I0/N,'r.');
legend('Two controls','Control U1 only','Without controls');
title('(D) I^*/N');
xlabel('Time t'); 
ylabel('I/N');

subplot(3,2,5)
plot(t,L2/N,'g--*',t,L21/N,'b-',t,L20/N,'r.');
legend('Two controls','Control U1 only','Without controls');
title('(E) L_2^*/N');
xlabel('Time t'); 
ylabel('L_2/N');

subplot(3,2,6)
plot(t,R/N,'g--*',t,R1/N,'b-',t,R0/N,'r.');
legend('Two controls','Control U1 only','Without controls', 'Location','southeast');
title('(F) R^*/N');
xlabel('Time t'); 
ylabel('R/N');


