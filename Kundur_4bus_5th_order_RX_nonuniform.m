
% R/X ratio
%ratio = [a b c];

% Drhoop Gains (%)
%K_p = 1*[1, 1, 1, 1];
%K_q = 1*[d, e, f, g];

% the number of inverters
n = 4;
% Base Peak Phase Voltage (V)
U_b = 230;
% Base Inverter Apparent Power
S_b = 10e3/3;
% Nominal Frequency (rad)
omega_0 = 2*pi*50;
% Filter Cut off Frequency (rad)
omega_c = 10*pi;
% Line Impedance (Ohm/km)
Z = 0.1587*(ratio+1i);
% Z_base (Ohm)
Z_b = U_b^2/S_b;
% Line Impedance (pu/km)
Z = Z/Z_b;

% Droop Gains
k_p = omega_0*diag(K_p./100);
k_q = diag(K_q./100);

% Kundur system line distances(km)
l = [6, 100, 3];
Z_set = Z.*l;

Y_n = [1/Z_set(1) -1/Z_set(1) 0 0;
    0 1/Z_set(2) -1/Z_set(2) 0
    0 0 1/Z_set(3) -1/Z_set(3)];

% Susceptance and Conductance matrices
B = -imag(Y_n);
G = real(Y_n);

% Incidence matrix
del = [1 0 0;
    -1 1 0;
    0 -1 1;
    0 0 -1];

% The state space matrix for Electro-Magnetic model
% s*x = A_5th*x, where x = (theta, omega, V, I_d, I_q)
A_5th = zeros(3*n+2*3);
% set of indixes for five state variables x = (theta, omega, V, I_d, I_q)
theta = 1:n;
omega = n+1:2*n;
V = 2*n+1:3*n;
I_d = 3*n+1:3*n+3;
I_q = 4*n+1:4*n+3;
% The state matrix construction
% trivial equation for s*theta = omega
A_5th(theta, omega) = eye(n);
% equation of the Low Pass filter s*omega = omega_c*omega
A_5th(omega, omega) = -omega_c*eye(n);
A_5th(omega, I_d) = -omega_c*k_p*del;
% equation of the Low Pass filter s*V = omega_c*V
A_5th(V, V) = -omega_c*eye(n);
A_5th(V, I_q) = omega_c*k_q*del;
A_5th(I_d, V) = omega_0*(diag(ratio)^2+eye(3))*B;
A_5th(I_d, I_d) = -omega_0*diag(ratio);
A_5th(I_d, I_q) = omega_0*eye(3);
A_5th(I_q, theta) = omega_0*(diag(ratio)^2+eye(3))*B;
A_5th(I_q, I_d) = -omega_0*eye(3);
A_5th(I_q, I_q) = -omega_0*diag(ratio);
