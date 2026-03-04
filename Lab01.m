%% LAB 01
format compact


%% Problem 1
clear all; close all; clc

% matrices + init cond.
A = [5 8; 1 3]
B = [4; -1]
C = [3 -4]
D = 7
x0 = [3; -3]

s = tf('s')

% Calculate inversion matrix
y = minreal(inv(s*eye(2)-A), 1e-3)

% Part a: ZERO-INPUT RESPONSE
X_zi = minreal(y * x0, 1e-3)
Y_zi = minreal(C * X_zi, 1e-3)

% Extract partial fractions to build analytical equations
% Residues - State 1: x1(t)
[num_x1a, den_x1a] = tfdata(X_zi(1), 'v')
[r_x1a, p_x1a] = residue(num_x1a, den_x1a)

% Residues - State 2: x2(t)
[num_x2a, den_x2a] = tfdata(X_zi(2), 'v')
[r_x2a, p_x2a] = residue(num_x2a, den_x2a)

% Residues - Output: y(t)
[num_ya, den_ya] = tfdata(Y_zi, 'v')
[r_ya, p_ya] = residue(num_ya, den_ya)

% Part b: Total Response to a step input
U = 4/s % laplace transform of a step input with amplitude 4
X_zs = minreal(y * B * U, 1e-3)
X_tot = minreal(X_zi + X_zs, 1e-3) % total zero-state reponse is X_zs(s) = y(s) * B * U(s)
Y_tot = minreal(C * X_tot + D * U, 1e-3)

% Extract partial fractions for part b
% Residues for total state 1: x1(t)
[num_x1b, den_x1b] = tfdata(X_tot(1), 'v')
[r_x1b, p_x1b] = residue(num_x1b, den_x1b)

% Residues for total state 2: x2(t)
[num_x2b, den_x2b] = tfdata(X_tot(2), 'v')
[r_x1b, p_x1b] = residue(num_x2b, den_x2b)

% Residues for output: y(t)
[num_yb, den_yb] = tfdata(Y_tot, 'v')
[r_yb, p_yb] = residue(num_yb, den_yb)


%% Problem 2
clear all; close all; clc

% matrices + init cond.
A = [0 1; -1 -1]
B= [27 0; -23 1]
C = [1 0]
x0 = [0; 0]

s = tf('s')
U = [0; 1] % input

% Q = inverse of (sI-A) 
Q = inv(s*eye(2)-A)
X = minreal(Q * (B*U + x0), 1e-3) % Total state response
Y = minreal(C*X, 1e-3) % Total output response

% Extract Partial Fractions and Euler Components
% Formula: 2*|R| * exp(Re(p)*t) * cos(Im(p)*t + angle(R))

% State 1: x1(t)
[num_x1, den_x1] = tfdata(X(1), 'v')
[r_x1, p_x1] = residue(num_x1, den_x1)

% Extract amplitude (2*|R|), phase (angle), and pole components
amp_x1 = 2 * abs(r_x1(1))      
phase_x1 = angle(r_x1(1))  
pole_re_x1 = real(p_x1(1))     
pole_im_x1 = imag(p_x1(1))     

% State 2: x2(t)
[num_x2, den_x2] = tfdata(X(2), 'v')
[r_x2, p_x2] = residue(num_x2, den_x2)

amp_x2 = 2 * abs(r_x2(1))      
phase_x2 = angle(r_x2(1))  
pole_re_x2 = real(p_x2(1))     
pole_im_x2 = imag(p_x2(1))   

% Output: y(t) 
% Since C = [1 0] and D = [0 0], y(t) mathematically equals x1(t).
[num_y, den_y] = tfdata(Y, 'v');
[r_y, p_y] = residue(num_y, den_y);

amp_y = 2 * abs(r_y(1))
phase_y = angle(r_y(1))


%% Problem 3
close all; clear all; clc

% matrices
A = [0 1 0 0; -1 -1 0 0; 0 0 -10 0; 0 0 0 -5];
B = [1; 0; 0; 1];
C = [1 0 1 1];
s = tf('s');
Q = minreal(inv(s*eye(4) - A), 1e-3);


% Part a: output response
x0_a = [0; 0; 0; 0];
U_a = (s + 1) / s^2;

% Calculate X&Y
X_a = minreal(Q * (B * U_a + x0_a), 1e-3);
Y_a = minreal(C * X_a, 1e-3) 

% extract poles and zeros
[num_ya, den_ya] = tfdata(Y_a, 'v');
[r_ya, p_ya] = residue(num_ya, den_ya)


% Part b: State response
x0_b = [0; 0; 1; 1];
U_b = 0;
X_b = minreal(Q * (B * U_b + x0_b), 1e-3)

% Extract and print r and p for state 3
[num_x3, den_x3] = tfdata(X_b(3), 'v');
[r_x3, p_x3] = residue(num_x3, den_x3)

% Extract and print r and p for state 4
[num_x4, den_x4] = tfdata(X_b(4), 'v');
[r_x4, p_x4] = residue(num_x4, den_x4)



%% Problem 4:
% matrices
A = [0 1 0 0; -1 0 0 0; 3 0 0 0; 0 0 2 0];
B = [1; 0; 0; 0];
C = [1 0 0 0];
x0 = [0; 0; 0; 0];

s = tf('s');
Q = minreal(inv(s*eye(4) - A), 1e-3);

%  PART a: OUTPUT RESPONSE TO u(t) = sin(t*sqrt(2))
% The Laplace transform of sin(w0*t) is w0 / (s^2 + w0^2) 
U_a = sqrt(2) / (s^2 + 2);

X_a = minreal(Q * (B * U_a + x0), 1e-3);
Y_a = minreal(C * X_a, 1e-3)
[num_ya, den_ya] = tfdata(Y_a, 'v');
[r_ya, p_ya] = residue(num_ya, den_ya)


%  PART b: OUTPUT RESPONSE TO u(t) = sin(t)
% Here w0 = 1, so the Laplace transform is 1 / (s^2 + 1)[cite: 279].
U_b = 1 / (s^2 + 1);

X_b = minreal(Q * (B * U_b + x0), 1e-3);
Y_b = minreal(C * X_b, 1e-3)
[num_yb, den_yb] = tfdata(Y_b, 'v');
[r_yb, p_yb] = residue(num_yb, den_yb)