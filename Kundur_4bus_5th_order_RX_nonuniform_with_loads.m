
% R/X ratio
ratio = 1.3;

% % Randomize R/X ratio 
% ratio = 0.4 + 2.1*rand(3,1);

% Drhoop Gains (%)
% K_p = 1.1881*[1; 1; 1; 1];
% K_q = K_p./0.3;
% 
K_p = [1; 1; 1; 1];
K_q = [1; 1; 1; 1];

% % Randomize K_q with fixed K_p
% K_q = K_p./(0.3+4.7*rand(4,1));

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
Z = Z./Z_b;

% Droop Gains
k_p = omega_0*diag(K_p./100);
k_q = diag(K_q./100);

% Kundur system line distances(km)
l = [6; 100; 3];
Z_line = Z.*l;

% Load impedances (Ohm)
Z_load = [20+1j; 25+1j;20+4.75j;40+12.58j];
Z_load = conj(Z_load)*1j;
% % Randomize load impedances within 10%
% Z_load = real(Z_load).*(1-0.5+1.5*rand(4,1))+1j*imag(Z_load).*(1-0.5+1.5*rand(4,1));
% Load impedsances (p.u.)
Z_load = Z_load./Z_b;

% Combine all impedances into one set
Z_set = [Z_line;Z_load];

% Y_n = [1/Z_set(1) -1/Z_set(1) 0 0;
%     0 1/Z_set(2) -1/Z_set(2) 0
%     0 0 1/Z_set(3) -1/Z_set(3)];
% 
% % Susceptance and Conductance matrices
% B = -imag(Y_n);
% G = real(Y_n);

% Incidence matrix
del_nw = [1 0 0;
    -1 1 0;
    0 -1 1;
    0 0 -1]';

del_load = eye(4);

del = [del_nw; del_load];
[n,m] = size(del');

L_p = 100/omega_0*diag(1./K_p);
L_q = 100*diag(1./K_q);


% The state space matrix for Electro-Magnetic model
% s*x = A_5th*x, where x = (theta, omega, V, I_d, I_q)
A_5th = zeros(3*n+2*m);
% set of indixes for five state variables x = (theta, omega, V, I_d, I_q)
theta = 1:n;
omega = n+1:2*n;
V = 2*n+1:3*n;
I_d = 3*n+1:3*n+m;
I_q = 3*n+m+1:3*n+2*m;
% The state matrix construction
% trivial equation for s*theta = omega
A_5th(theta, omega) = eye(n);
% equation of the Low Pass filter s*omega = omega_c*omega
A_5th(omega, omega) = -omega_c*L_p;
A_5th(omega, I_d) = -omega_c*del';
% equation of the Low Pass filter s*V = omega_c*V
A_5th(V, V) = -omega_c*L_q;
A_5th(V, I_q) = omega_c*del';
A_5th(I_d, V) = omega_0*del;
A_5th(I_d, I_d) = -omega_0*diag(real(Z_set));
A_5th(I_d, I_q) = omega_0*diag(imag(Z_set));
A_5th(I_q, theta) = omega_0*del;
A_5th(I_q, I_d) = -omega_0*diag(imag(Z_set));
A_5th(I_q, I_q) = -omega_0*diag(real(Z_set));
Gamma = eye(3*n+2*m);
Gamma(omega,omega) = L_p;
Gamma(V,V) = L_q;
Gamma(I_d,I_d) = diag(imag(Z_set));
Gamma(I_q,I_q) = diag(imag(Z_set));

% Caclculate spectrum
[W,L] = eig(A_5th,Gamma);
L = diag(L);
W = W(:,abs(L)<1e10);
L = L(abs(L)<1e10);
[~,ind] = sort(real(L));
L = L(ind);
W = W(:,ind);
% Dominant oscillatory modes
W_dom = W(:, imag(L)>1e-10&abs(L)<500);
L_dom = L(imag(L)>1e-10&abs(L)<500);
