
% Extended Kalman Filter
function [x, p, y, S] = EKF(...
    n, y, Q, R_meas, mu0, P0, Fx, Omega, Hx, Dt, rho, drho, phi, mu)

    [p, nS, T] = size(y);
    
    x_t = zeros(n, T);
    P_t = zeros(n, n, T);
    
    y_m = NaN * zeros(p, nS, T);
    
    
    mu_prev = mu0;
    P_prev = P0;
    
    x_t(:, 1) = mu_prev;
    P_t(:, :, 1) = P_prev;
    
    k = (0:(T-1))';
    
    tspan = [0, 1, Dt];
    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    
    S_t = NaN * zeros(p*nS, p*nS, T);
    
    % Find first S
    yinds = find(~isnan(y(1, :, 1)));
    n_meas = length(yinds);
    if n_meas > 0 
        H = zeros(n_meas*p, n);
        for j = 1:n_meas
            H((p*j-2):(p*j),:) = Hx(k(1), mu_prev, yinds(j));
        end
        R = kron(eye(n_meas), R_meas);
        S_t(1:n_meas*p, 1:n_meas*p, 1) = H*P0*H'+R;
    end

    for i = 2:T

        [~, x] = ode45(@(t, S) nl_sim_nonoise(t, S, mu),...
            tspan, mu_prev, opts);
        
        x = x(end, :);
        
        F = Fx(mu_prev);
        P = F*P_prev*F' + Omega*Q*Omega';
        
        yframe = y(:, :, i);
        yinds = find(~isnan(yframe(1, :)));
        n_meas = length(yinds);
        mu_prev = x';
        P_prev = P;
        if n_meas > 0
            R = kron(eye(n_meas), R_meas);
            H = zeros(p*n_meas, n);
            for j = 1:n_meas
                H((p*j-2):(p*j), :) = Hx(k(i), x, yinds(j)); 
            end
            
            yi = yframe(:, yinds);
            yi = yi(:);
            
            rho_mi = rho(k(i)*Dt, x, yinds);
            drho_mi = drho(k(i)*Dt, x, yinds);
            phi_mi = phi(k(i)*Dt, x, yinds);
            ymi2 = [rho_mi; drho_mi; phi_mi];
            ymi = ymi2(:);
            
            e_yi = yi - ymi;
            
            S = H*P*H'+R;
            K = P*H'*(S\eye(p*n_meas));
            mu_prev = x' + K*e_yi;
            P_prev = (eye(n)-K*H)*P;
            
            S_t(1:p*n_meas, 1:p*n_meas, i) = S;
            
            y_m(:, yinds, i) = ymi2;
        end
        
        x_t(:, i) = mu_prev;
        P_t(:, :, i) = P_prev;
        
    end

end