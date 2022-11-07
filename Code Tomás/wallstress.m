% Wall stress
function tau = wallstress(rho, V, rho_d, V_d, T_d, r_d)
% n = 0.925e-5*(T_d./300).^1.1;
% Cd  =0.002;
% tau_d = 12.*n.*V_d./(2.*r_d) + 2*Cd.*rho_d.*V_d.^2;
% tau = tau_d./(rho_ref.*a0.^2); % non-dimensional

%f = 0.00028; % average value from Gasdynamics (needs improvement)
mu = 0.925*10^-5.*(T_d./300).^1.1;
Re = rho_d.*V_d.*2.*r_d./mu;
f = (0.79*log(Re)-1.64).^-2;
tau = 0.5.*f.*rho.*V.^2;%[non-dimensional]
%tau = tau_d./(rho_ref.*a0.^2); % non-dimensional

end