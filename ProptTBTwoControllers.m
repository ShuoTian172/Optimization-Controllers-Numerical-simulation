toms t
p = tomPhase('p', t, 0, 5, 50);
setPhase(p);
tomStates S L1 I L2 R
tomControls U1 U2

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


% Initial guess
% Note: The guess for t_f must appear in the list before expression involving t.
x0 = {collocate(U2==0.001),...
collocate(U1==0.001)};

% Box constraints
cbox = {0 <= collocate(U1) <= 1
0 <= collocate(U2) <= 1};

% Boundary constraints
cbnd = {initial({S == N*76/120; L1 == N*36/120
I == N*5/120; L2 == N*5/120; R == N*1/120})};

% ODEs and path constraints
ceq = collocate({
dot(S) == miu*N - b/N*I*S - miu*S;
dot(L1) == b/N*I*(S + rho*L2 + rhor*R) - (d + tao1 + miu)* L1;
dot(I) == phi*d*L1 + w*L2 + wr*R - (tao0 + e1*U1 + miu)*I;
dot(L2) == (1 - phi)*d*L1 - rho*b/N*I*L2 - (w + e2*U2 + tao2 + miu)*L2;
dot(R) == (tao0 + e1*U1)*I + tao1*L1 + (tao2 + e2*U2)*L2 - rhor*b/N*I*R-(wr + miu)*R});


% Objective
objective = integrate(I + L2 + 0.5*W1*U1^2 + 0.5*W2*U2^2);

% Solve the problem
options = struct;
solution = ezsolve(objective, {cbox, cbnd, ceq}, x0, options);
t = subs(collocate(t),solution);
S = subs(collocate(S),solution);
L1 = subs(collocate(L1),solution);
I = subs(collocate(I),solution);
L2 = subs(collocate(L2),solution);
R = subs(collocate(R),solution);
U1 = subs(collocate(U1),solution);
U2 = subs(collocate(U2),solution);

subplot(3,2,1)
plot(t,(I+L2)/N,'b-',t,(I1+L21)/N,'g*-',t,(I0+L20)/N);
legend('Two controls','Control U1 only','Without controls');
title('(A) (I^*+L_2^*)/N');
xlabel('Time t'); 
ylabel('(I+L2)/N');

subplot(3,2,2)
plot(t,S/N,'b-',t,S1/N,'g--*',t,S0/N);
legend('Two controls','Control U1 only','Without controls');
title('(B) S^*/N');
xlabel('Time t'); 
ylabel('S/N');

subplot(3,2,3)
plot(t,L1/N,'b-',t,L11/N,'g*-',t,L10/N);
legend('Two controls','Control U1 only','Without controls');
title('(C) L_1^*/N');
xlabel('Time t'); 
ylabel('L_1/N');

subplot(3,2,4)
plot(t,I/N,'b-',t,I1/N,'g*-',t,I0/N);
legend('Two controls','Control U1 only','Without controls');
title('(D) I^*/N');
xlabel('Time t'); 
ylabel('I/N');

subplot(3,2,5)
plot(t,L2/N,'b-',t,L21/N,'g*-',t,L20/N);
legend('Two controls','Control U1 only','Without controls');
title('(E) L_2^*/N');
xlabel('Time t'); 
ylabel('L_2/N');

subplot(3,2,6)
plot(t,R/N,'b-',t,R1/N,'g*-',t,R0/N);
legend('Two controls','Control U1 only','Without controls');
title('(F) R^*/N');
xlabel('Time t'); 
ylabel('R/N');


