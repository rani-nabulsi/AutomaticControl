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