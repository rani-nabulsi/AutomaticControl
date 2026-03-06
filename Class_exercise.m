format compact

%% Problem 1:

% System Matrices
A = [0 1; -2 -3];
B = [1; 0];

% Initial Conditions + Laplace variable
x0 = [2; 2];
s = tf('s');

% Input Signal
U = 2/s;

% Compute invrse matrix
n = 2;
Q = inv(s*eye(n) - A);  %% inv computes the inverse
x = Q*(x0+B*U);
x = minreal(x, 1e-3);

%% PFE -> residue
[num, den] = tfdata(x(1), 'v');
[residues, poles] = residue(num, den) % PFE



%% Problem 2:
close all; clear all; clc
A = [-3 2; -2 -3]
B = [1; 0]
C = [0 1]
D = 0

x0 = [1; 1];
s = tf('s');
U = 1/s;
Q = inv(s*eye(2)-A);
X = minreal(Q * (x0 + B*U), 1e-3)
Y = X(2)
[num,den] = tfdata(Y,'v')
[r,p] = residue(num,den)

M = abs(r(1)) % modulus |R|
2*M %2|R|

phi = angle(r(1)) % phase


%% Problem 3: Transfer Function
clear all; close all; clc
A = [-3 2; -2 -3]; B = [1; 0]; C = [0 1]; D = 0;
s = tf('s')

Q = inv(s*eye(2)-A)
H = minreal(C*Q*B+D, 1e-3) % matrix inversion with polynomials
% H = tf(ss(A,B,C,D)) % alternative, can be defined as a transfer given is its ss

% roots and poles
zero(H) 
pole(H) 



%% Problem 4: Output Response of LTI using Laplace Transforms
clear all; close all; clc
s = tf('s')
H = minreal((2*s+1)/(s+4)^2, 1e-3) % matrix inversion (wrap to cancel out overlapping)
U = 2/s^2 % input

% Calculate Output Y(s)
Y = minreal(H*U, 1e-3)

% Extract Numerator & Denominator Polynomials
[num,den] = tfdata(Y, 'v');
[r, p] = residue(num, den)


%{
    Theorem 1
An LTI System is internally stable if and only if all its natural
modes are bounded.

    Theorem 2
An LTI system is asymptotically stable if and only if all its natural
modes are convergent.

    Theorem 3
An LTI system is internally unstable if and only if there is at least one
divergent natural mode.

    Theorem 4
An LTI system is internally stable if and only if all the eigenvalue of
A have real part <= 0 and those with real part = 0 have minimal
multiplicity

    Theorem 5
An LTI system is asymptotically stable if and only if all the eigenvalues
of A have real part < 0

    Theorem 6
An LTI system is internally unstable if there exists either an eigenvalue
with real part > 0 or an eignvalue with real part = 0 and minimal
multiplicity
%}

%% Problem 5: stability of LTI systems
clear all; close all; clc
A = [0 1 0; 0 0 0; 0 0 -1]
eig(A)
s = tf('s')

% METHOD: compute (sI-A)^-1
% Calculate resolvent matrix
% zpk(...): formats the output into Zero-Pole-Gain format
zpk(minreal(inv(s*eye(3)-A), 1e-3))



%% Problem 6
clear all; close all; clc
A = [0 1 0; 0 0 0; 0 0 -1]
eig(A)
s = tf('s')
zpk(minreal(inv(s*eye(3)-A), 1e-3))
roots(minpoly(A))