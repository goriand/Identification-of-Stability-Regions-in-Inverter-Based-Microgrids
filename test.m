for k1=0.3:0.1:5
    k=500;
    step=k;
    for i = 1:log2(500 / 1e-2)
        step=step/2;
        k_p = k*ones(1,10);
        k_q = k/k1*ones(1,10);
        K_p(inv_ind) = k_p;
        K_q(inv_ind) = k_q;
        L_p = 100/omega_0*diag(1./K_p);
        L_q = 100*diag(1./K_q);
        % The state matrix update
        % equation of the Low Pass filter s*omega = omega_c*omega
        A_5th(omega, omega) = -omega_c*L_p;
        % equation of the Low Pass filter s*V = omega_c*V
        A_5th(V, V) = -omega_c*L_q;
        Gamma(omega,omega) = L_p;
        Gamma(V,V) = L_q;
        L = eig(A_5th,Gamma);
        L = L(abs(L)<1e10);
%         % Dominant oscillatory modes
%         L_dom = L(imag(L)>1e-10&abs(L)<500);
%         IEEE123bus_nonuniform_rx
        %Kundur_4bus_5th_order_RX_nonuniform_with_loads
        if all(real(L) < 0)
            k=k+step;
        else
            k=k-step;
        end
    end
    mu_percent102_k1_k5(int64(k1/0.1 - 2)) = k;
end

