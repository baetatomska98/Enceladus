function [dR_dx, rho_eq] = dR_dx_MC(T_d,rho_d,M,b)
% Computation of the growth rate of the icy grains, see Schmidt et al
%Rg = 461.52; %Gas constant water vapour [J/(kg K)]
rho_grain = 920; %density of ice [kg/m3]
g = 4/3;
m0 = 2.992*10^(-26);
kb = 1.38065*10^(-23);
%rho_eq = p_eq_sg(T)./(Rg*T);
rho_eq = p_eq_sg(T_d).*m0./(kb*T_d);
%dRdx = beta./(sqrt(2*pi*g).*rho_grain).*(rho - rho_eq).*(rho > rho_eq)./M;

dR_dx = b./(sqrt(2*pi*g).*rho_grain).*(rho_d - rho_eq)./M;
% If rho < rho_eq -> EVAPORATION
% Energy release?