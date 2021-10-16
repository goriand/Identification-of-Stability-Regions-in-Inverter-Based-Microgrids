% % Initialize the network parameters
%IEEE_123_config_lines_initialization

% remove isolated nodes
nabla( :, ~any(nabla,1) ) = [];

% incidence matrix of the network graph
del = nabla;
% n - the number of inverters
% m - the number of buses
[m,n] = size(del);

% Inverter locations
%inv_ind = [19 51 65 89 103];
inv_ind = [95   115    79     5   102   112    81    91    89    47];
%inv_ind = [95   36   38];
%inv_ind = [89  91];

% % % Drhoop Gains (%)
% % k_p = k*ones(size(inv_ind));
% % k_q = k/k1*ones(size(inv_ind));
% 
% k_p = 1*ones(size(inv_ind));
% 
% % %enlarge stability regions using Gershgorin's disks
% % k_p = [26.7548   20.5420   26.1948   21.1612   25.5604   49.8624   32.9895    4.8764    5.3101 67.0016];
k_p = 100*[0.2650
    0.2034
    0.2594
    0.2096
    0.2532
    0.4938
    0.3267
    0.0483
    0.0526
    0.6636]';
% % 
% % % Uniform droops
% k_p = 5.29*ones(size(inv_ind));
% 
% k_q = 1/0.3*k_p;
% 
% % %random k_q with fixed k_p
% % k_q = k_p.*(1/5 + (1/0.3 - 1/5)*rand(size(inv_ind)));
k_p = 10*ones(1,10);
k_q = 4*ones(1,10);

% Sensitivity calculation
% Change of droop gain by epsilon
% k_p(o) = k_p(o) + epsilon;
% k_q(o) = k_q(o) + k/k1*epsilon;

K_p = +Inf*ones(1,n);
K_q = +Inf*ones(1,n);
K_p(inv_ind) = k_p;
K_q(inv_ind) = k_q;

% Base Peak Phase Voltage (V)
U_b = 120; 
% Base Inverter Apparent Power
S_b = 5*10e3/3;
% Nominal Frequency (rad)
omega_0 = 2*pi*50;
% Filter Cut off Frequency (rad)
omega_c = 10*pi;


% Z_base (Ohm)
Z_b = U_b^2/S_b;


% Inverse Droop Gains
L_p = 100/omega_0*diag(1./K_p);
L_q = 100*diag(1./K_q);


% Load Impedances (pu)
Z_load = S_b./(Pload_balanced-1j*Qload_balanced);

%% Implement random R with fixed X for the lines
% Uniform distribution of R in [0.4:4] of lines
% a = 0.4;
% b=4;
% Z_line = (a + (b-a)*rand(size(Z_line)) + 1j).*imag(Z_line);
% 
% % Uniform distribution of R in loads
% Z_load = (0.01+10*rand(size(Z_load))+1j).*imag(Z_load);

% Line Impedance (pu)
Z_set = [Z_line /Z_b; Z_load];

% % Implement homogeneous R/X ratio
% ratio = 1.3;
% % Fix line reactances X
% Z_set = (ratio+1j)*imag(Z_set);

%% The state space matrix for Electro-Magnetic model
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

%% Caclculate spectrum
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

% %% Addmitance matrix
% Ybus = nabla'*diag(1./Z_set)*nabla;
% % Kron reduction
% noinv_ind = setdiff(1:n,inv_ind);
% Ybus_red = Ybus(inv_ind,inv_ind) - ...
%     Ybus(inv_ind,noinv_ind)/Ybus(noinv_ind,noinv_ind)*Ybus(noinv_ind,inv_ind);
% 
% % 1/X Laplacian matrix
% B = nabla'*diag(1./imag(Z_set))*nabla;
% B_red = B(inv_ind,inv_ind) - ...
% B(inv_ind,noinv_ind)/B(noinv_ind,noinv_ind)*B(noinv_ind,inv_ind);
% % B' 'Droop-Weighted' Laplacian matrix
% B_dash = diag(K_p(inv_ind)./100)*B_red;
% [PSI,M, PHI] = eig(B_dash);
% [~,ind] = sort(diag(M));
% M = diag(M(ind,ind));
% PSI = PSI(:,ind);
% PHI = PHI(:,ind);
% 
% % Eigenvalue mu sensitivities
% for i=1:length(inv_ind)
%     dmu(i) = PSI(i,end)*B_red(i,:)*PSI(:,end)/(PSI(:,end)'*PHI(:,end))/100;
% end