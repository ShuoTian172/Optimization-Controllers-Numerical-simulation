toms t
p = tomPhase('p', t, 0, 5, 50);
setPhase(p);
tomStates S1 L11 I1 L21 R1
tomControls U11

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

U21=0;

W1 = 500;
W2 = 50;


% Initial guess
% Note: The guess for t_f must appear in the list before expression involving t.
x0 = {collocate(U11==0.001)};

% Box constraints
cbox = {0 <= collocate(U11) <= 1};

% Boundary constraints
cbnd = {initial({S1 == N*76/120; L11 == N*36/120
I1 == N*5/120; L21 == N*5/120; R1 == N*1/120})};

% ODEs and path constraints
ceq = collocate({
dot(S1) == miu*N - b/N*I1*S1 - miu*S1;
dot(L11) == b/N*I1*(S1 + rho*L21 + rhor*R1) - (d + tao1 + miu)* L11;
dot(I1) == phi*d*L11 + w*L21 + wr*R1 - (tao0 + e1*U11 + miu)*I1;
dot(L21) == (1 - phi)*d*L11 - rho*b/N*I1*L21 - (w + e2*U21 + tao2 + miu)*L21;
dot(R1) == (tao0 + e1*U11)*I1 + tao1*L11 + (tao2 + e2*U21)*L21 - rhor*b/N*I1*R1-(wr + miu)*R1});


% Objective
objective = integrate(I1 + L21 + 0.5*W1*U11^2);

% Solve the problem
options = struct;
solution = ezsolve(objective, {cbox, cbnd, ceq}, x0, options);
t = subs(collocate(t),solution);
S1 = subs(collocate(S1),solution);
L11 = subs(collocate(L11),solution);
I1 = subs(collocate(I1),solution);
L21 = subs(collocate(L21),solution);
R1 = subs(collocate(R1),solution);
U11 = subs(collocate(U11),solution);
% U21 = subs(collocate(U21),solution);



subplot(3,2,1)
plot(t,(I+L2)/N,'g--*',t,(I1+L21)/N,'b-');
legend('Two controls','Control U1 only');
title('(A) (I^*+L_2^*)/N');
xlabel('Time t'); 
ylabel('(I+L2)/N');

subplot(3,2,2)
plot(t,S/N,'g--*',t,S1/N,'b-');
legend('Two controls','Control U1 only');
title('(B) S^*/N');
xlabel('Time t'); 
ylabel('S/N');

subplot(3,2,3)
plot(t,L1/N,'g--*',t,L11/N,'b-');
legend('Two controls','Control U1 only');
title('(C) L_1^*/N');
xlabel('Time t'); 
ylabel('L_1/N');

subplot(3,2,4)
plot(t,I/N,'g--*',t,I1/N,'b-');
legend('Two controls','Control U1 only');
title('(D) I^*/N');
xlabel('Time t'); 
ylabel('I/N');

subplot(3,2,5)
plot(t,L2/N,'g--*',t,L21/N,'b-');
legend('Two controls','Control U1 only');
title('(E) L_2^*/N');
xlabel('Time t'); 
ylabel('L_2/N');

subplot(3,2,6)
plot(t,R/N,'g--*',t,R1/N,'b-');
legend('Two controls','Control U1 only');
title('(F) R^*/N');
xlabel('Time t'); 
ylabel('R/N');

