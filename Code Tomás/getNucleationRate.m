% getNucleationRate %
function J_H2O = getNucleationRate(T, rho, Rg)
% Equations are from Wölk & Strey 2001 or 2003 paper.
water_molarmass = 18.01528; % g/mol 
N_A = 6.022140857*10^23; % Avogadro constant, molecules/mol
m_m = water_molarmass/N_A; % gram per molecule
m_m = m_m*10^-3; %kg/molecule
k = 1.38064852*10^(-23); % Boltzmann constant [m^2kgs^-2K^-1]
T_C = 647.15; % 647.15? Cite: n.4 (DataPaper_WölkStrey)
x = (T-225)./46.2; % DataPaper
t_r = (T_C - T)./T_C; % DataPaper
rho_vm = 0.08*tanh(x) + 0.7415*t_r.^0.33 + 0.32; % Density [g/cm3] for v_m calc.
rho_vm = rho_vm*10^3;
v_m = m_m./rho_vm; % Molecular volume of water [1/m^3]
sigma = (93.6635 + 0.009133*T - 0.000275*T.^2)*10^(-3); % cite: n.4,7 (DataPaper_WölkStray)
p_eq = exp(77.34491-7235.42465./T - 8.2.*log(T) + 0.0057113.*T); % cite: n.66 (DataPaper_WölkStrey)
rho_eq = p_eq_lg(T)./(Rg*T);
S = rho./rho_eq;
p_v = p_eq .* S;
J_DB = sqrt((2.*sigma)./(pi.*m_m)).*v_m.*(p_v ./ (k.*T)).^2.*exp((-16.*pi.*v_m.^2.*sigma.^3)./(3.*(k.*T).^3 .* log(S).^2)); 
J_H2O = J_DB .* exp(-27.56 + 6.5*10^3./T); %[1/m^3s]
J_H2O = J_H2O./(10^6); % %[1/cm^3s]

end