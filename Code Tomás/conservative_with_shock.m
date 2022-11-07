function [M, V_d, p_d, m_d, x_d, rho_d, T_d] = conservative_with_shock(x, dx, n, C, A, nt, kg, throat, Cx, Rg, pe, L, cv, T_ref, rho_ref, p_ref, r, c)
global cp crackLength circShape gn E Friction U1 U2 U3 b simp M_ref convergenceCriteria critRes rho_m p_m T_m M_m f R fp As A_d r_d dR_dx Edot_L As A_d r_d stop
rho_ice = 920;                      % Ice density [kg/m3]
H2O_mol_mass = 0.01801528;          % Molar mass H2O [kg/mol]
NA = 6.02214076*10^(23);            % nr. of molecules in one mole
L_h = 2.836e6;                      % Latent Heat [J/kg]
% Throat location (dimensional)
x_d = x .* L;
% from Mini_moon.m (non-dimensional) %
rho = rho_m;
T = T_m;
p = p_m;
M = M_m;
V = M .* sqrt(T);

% Solution vectors
U1 = rho .* A;
U2 = rho .* A .* V;
U3 = rho .* A .* (T / (kg - 1) + kg / 2 * V.^2);

for k = 1 : nt    
    % Copies of solution vectors
    U1_old = U1;
    U2_old = U2;
    U3_old = U3;

    % dt calculations
    dt = min(C * dx./(sqrt(T) + V));

    % Flux vectors
    F1 = U2;
    F2 = U2.^2./U1 + (kg - 1)/kg * (U3 - kg * U2.^2./(2 * U1));
    F3 = kg * U2 .* U3./U1 - kg * (kg - 1)/2 * U2.^3./U1.^2;
    [rho_d, T_d, p_d, V_d, A_d, x_d, dx_d, r_d] = dimensionalize(rho, T, p, M, kg, Rg, x, dx, rho_ref, T_ref, p_ref, L, r, As, circShape, crackLength);
    Q = mean(rho_d .* A_d .* V_d); % [kg/s], average because; mass flow = C
    dR_dx = dR_dx_MC(T_d, rho_d, M, b); % [-]
    R = getParticleRadius(dR_dx, x_d); % [m], [-]
    fp = solid_frac_MC(Q, dx_d, T_d, rho_d, R, dR_dx, A_d, Rg); % uses rho_d(i), T_d(i), M(i)
    % Updating solid fraction, Change to Simpson's rule  %
    f = getSolidFraction(fp, x_d);
    % Non-dimensionalize
    [fp_nd, Lh_nd, Rg_nd, E_nd] = nondimensionalize(fp, L, L_h, cv, T_ref, Rg, E, rho_ref, kg);
    % Calculate energy of latent heat in all segments
    Edot_L = Edot_source(rho, f , fp, dx_d, V, L_h, A, cv, Rg, dx, fp_nd, Lh_nd, x, T, kg); % Latent heat as Source term
    mdot_d  = mdot_sink(f, rho, A, V); % Mass reduction ~Phase change vapor -> ice
    if Friction == true
       tau = wallstress(rho, V, rho_d, V_d, T_d, r_d); % [-] [kg/ms2] 
       friction = getfriction(tau, c, dx); % Non-dimensional
    end
    % Predictor step:
    % Compute derivatives of the solution vectors
    for i = 2 : n - 1
        J2(i) = 1 / kg * rho(i) .* T(i) * fd_dx(A, dx, i);
        dU1dt_p(i) = -fd_dx(F1, dx, i) - fd_dx(mdot_d, dx, i); % mass is reduced ~Deposition
        if Friction == true
            dU2dt_p(i) = -fd_dx(F2, dx, i) + J2(i) - V(i) * rd_dx(mdot_d, dx, i) - friction(i)/dx;
        else
            dU2dt_p(i) = -fd_dx(F2, dx, i) + J2(i) - V(i) * rd_dx(mdot_d, dx, i);
        end
        dU3dt_p(i) = -fd_dx(F3, dx, i) + fd_dx(Edot_L, dx, i);
        latentheat_fd(i) = fd_dx(Edot_L, dx, i);
        % Shock Capturing equations, based on 'old' (current) U1-3
        cc = Cx * abs(cdiff(p,i))/(p(i+1) + 2 * p(i) + p(i-1));
        S1(i) = cc * cdiff(U1,i);
        S2(i) = cc * cdiff(U2,i);
        S3(i) = cc * cdiff(U3,i);

        % PREDICTED solution vectors (including S-terms)
        U1(i)= U1(i) + dU1dt_p(i) * dt + S1(i);
        U2(i)= U2(i) + dU2dt_p(i) * dt + S2(i);
        U3(i)= U3(i) + dU3dt_p(i) * dt + S3(i);
    end
    % PREDICTED primitive variables and flux vectors
    rho = U1 ./ A;
    T = (kg - 1) * (U3 ./ U1 - kg / 2 * (U2 ./ U1).^2);
    p = rho .* T;
    % Compute new M %
    M = V ./ sqrt(T);
    [rho_d, T_d, p_d, V_d, A_d, x_d, dx_d, r_d] = dimensionalize(rho, T, p, M, kg, Rg, x, dx, rho_ref, T_ref, p_ref, L, r, As, circShape, crackLength);

    % Dimensionalize for computation: f', dE, dR. ('_d' = dimensionalized) %
    Q = mean(rho_d .* A_d .* V_d); % average because; mass flow = C
    dR_dx = dR_dx_MC(T_d,rho_d,M, b);
    R = getParticleRadius(dR_dx, x_d);
    [fp, gn] = solid_frac_MC(Q, dx_d, T_d, rho_d, R, dR_dx, A_d, Rg); % uses rho_d(i), T_d(i), M(i)
    % Updating solid fraction, Change to Simpson's rule %
    f = getSolidFraction(fp, x_d);
    % Non-dimensionalize
    [fp_nd, Lh_nd, Rg_nd, E_nd] = nondimensionalize(fp, L, L_h, cv, T_ref, Rg, E, rho_ref, kg);
    % Calculate energy of latent heat in all segments
    Edot_L = Edot_source(rho, f , fp, dx_d, V, L_h, A, cv, Rg, dx, fp_nd, Lh_nd, x, T, kg); % Latent heat as Source term
    mdot_d  = mdot_sink(f, rho, A, V); % Mass reduction ~Deposition
    if Friction == true
       tau = wallstress(rho, V, rho_d, V_d, T_d, r_d); % [-] [kg/ms2] 
       friction = getfriction(tau, c, dx); % Non-dimensional
    end
    F1 = U2;
    F2 = U2.^2 ./ U1 + (kg - 1) / kg * (U3 - kg * U2.^2 ./ (2 * U1));
    F3 = kg * U2 .* U3 ./ U1 - kg * (kg - 1) / 2 * U2.^3 ./ U1.^2;
    % Corrector step
    for j = 2 : n - 1
        J2(j) = 1 / kg * rho(j) .* T(j) * rd_dx(A, dx, j);
        dU1dt_c(j) = -rd_dx(F1, dx, j) - rd_dx(mdot_d, dx, j); % mass is reduced ~Deposition
        if Friction == true
            dU2dt_c(i) = -fd_dx(F2, dx, i) + J2(i) - V(i) * rd_dx(mdot_d, dx, i) - friction(i)/dx;
        else
            dU2dt_c(i) = -fd_dx(F2, dx, i) + J2(i) - V(i) * rd_dx(mdot_d, dx, i);
        end
        dU3dt_c(j) = -rd_dx(F3, dx, j) + rd_dx(Edot_L, dx, j);
        latentheat_rd(j) = rd_dx(Edot_L, dx, j);

        % Shock Capturing equations
        cc = Cx * abs(cdiff(p,j))/(p(j+1) + 2 * p(j) + p(j-1));
        S1(j) = cc * cdiff(U1,j);
        S2(j) = cc * cdiff(U2,j);
        S3(j) = cc * cdiff(U3,j);
    end
    % Averaging corrector step and predictor step
    dU1dt_a = 0.5 * (dU1dt_p + dU1dt_c);
    dU2dt_a = 0.5 * (dU2dt_p + dU2dt_c);
    dU3dt_a = 0.5 * (dU3dt_p + dU3dt_c);

    % Updating solution vectors after averaging
    for m = 2 : n - 1
        U1(m) = U1_old(m) + dU1dt_a(m) * dt + S1(m);
        U2(m) = U2_old(m) + dU2dt_a(m) * dt + S2(m);
        U3(m) = U3_old(m) + dU3dt_a(m) * dt + S3(m);
    end   
    resU1(k) = sum(abs(U1 - U1_old))/sum(U1_old);
    resU2(k) = sum(abs(U2 - U2_old))/sum(U2_old);
    resU3(k) = sum(abs(U3 - U3_old))/sum(U3_old);

    % Boundary conditions
    % Inlet, fix rho & T
    U1(1) = rho(1)*A(1);
    U2(1) = 2*U2(2) - U2(3);
    V(1) = U2(1)/U1(1);
    T(1) = 1;
    U3(1) = U1(1) * (T(1) / (kg - 1) + kg / 2 * V(1)^2);

    % Outlet
    U1(n) = 2*U1(n-1) - U1(n-2);
    U2(n) = 2*U2(n-1) - U2(n-2);
    if pe > 0 % fix for shock-capturing
        U3(n) = pe * A(n)/(kg - 1) + kg / 2 * U2(n) * (U2(n) / U1(n));%V(n);
        p(n) = pe;
    else % if pe = undefined [subsonic-supersonic],
        U3(n) = 2*U3(n-1) - U3(n-2); %Only when pe/p0 is NOT set
    end
    % Updation of flow field variables
    rho = U1 ./ A;
    V = U2 ./ U1;
    T = (kg - 1) * (U3 ./ U1 - kg / 2 * (U2 ./ U1).^2);
    T(1) = 1;
    p = rho .* T;

    % Defining mass flow and mach number
    m = rho .* A .* V;
    m_d = rho_d .* A_d .* V_d;
    M = V ./ sqrt(T);
    % Dimensional for computation: f', dE, dR. ('_d' = dimensionalized) %
    [rho_d, T_d, p_d, V_d, A_d, x_d, dx_d, r_d] = dimensionalize(rho, T, p, M, kg, Rg, x, dx, rho_ref, T_ref, p_ref, L, r, As, circShape, crackLength);

    % Capturing the values at inlet, throat and exit for each iterations
    mach_in(k) = M(1);
    pressure_in(k) = p(1);
    density_in(k) = rho(1);
    temperature_in(k) = T(1);

    mach_thr(k) = M(throat);
    pressure_thr(k) = p(throat);
    density_thr(k) = rho(throat);
    temperature_thr(k) = T(throat);
    mach_ex(k) = M(n);
    pressure_ex(k) = p(n);
    density_ex(k) = rho(n);
    temperature_ex(k) = T(n);
    gn_max = max(gn);
    gn_max_n = find(gn_max == gn);
    gn_max_z(k) = x_d(gn_max_n);
    f_exit(k) = f(end);

    if convergenceCriteria == true && resU1(k) < critRes && resU2(k) < critRes  && resU3(k) < critRes
        break % Exit for-loop when convergence criteria is met
    end
    if isreal(M) == 0 || any(T<0) == true || sum(isnan(M)) > 0 || sum(isnan(T)) > 0
        display(['Imaginary numbers at ', num2str(k), ' iterations.'])
        break
    end
    if convergenceCriteria == true && resU1(k) < critRes && resU2(k) < critRes  && resU3(k) < critRes
        break % Exit for-loop when convergence criteria is met
    end
    if stop == true % type in command window "stop = true" to exit loop
        break
    end
end
[rho_d, T_d, p_d, V_d, A_d, x_d, ~, r_d] = dimensionalize(rho, T, p, M, kg, Rg, x, dx, rho_ref, T_ref, p_ref, L, r, As, circShape, crackLength);
end