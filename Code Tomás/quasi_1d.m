% Quasi 1D nozzle flows
clear all
close all
clc
global cp crackLength circShape gn E Friction U1 U2 U3 b simp M_ref convergenceCriteria critRes rho_m p_m T_m M_m f R fp As A_d r_d dR_dx Edot_L As A_d r_d stop
C = 0.001;                                          % Courant number, reduces dt and thus dudt. Set as low as possible, while not oscillating
Cx = 0.2;                                           % Viscosity constant (Adjustable)
nt = 70000;                                         % Iterations / Total time steps
n = 201;                                            % Number of segments
%%% Constants %%%
cv = 1384.5;                                        % specific heat capacity of water vapor at constant volume [J/(kg K)]
cp = 1846;                                          % specific heat capacity of water vapor at constant pressure [J/(kg K)]
m0 = 2.988e-26;                                     % Mass of a water molecule [kg]
kb = 1.38065e-23;                                   % Boltzmann constant [J/K]
Rg = 461.5;                                         % Specific gas constant [J/kgK]
%%% Reservoir conditions @Triple Point, from Schmidt %%%
pe = 0;
T_res = 273.16;
p_res = 611.2;
rho_res = p_res/(Rg*T_res);
kg = 4/3;                   % Specific heat ratio
b = 0.2;                    % Condensation coefficient
critRes = 10^(-10);         % Convergence criteria for the residuals
crackLength = 200;          % Crack length [m]
circShape = false;          % Circular shape channel or rectangular
Friction = false;
convergenceCriteria = true; % Terminate run at convergence criteria
simp = true;                % Integration with simpson method

L = 150;                     % Nozzle length [m]
x = linspace(0,L,n);
% r = zeros(1,length(x));
% x_th = 0.75;
% r_max = 0.07;
% r_th = 0.05;
% r(1:100) = r_max - (r_max - r_th)/x_th * x(1:100);
% r(101:end) = r_th + (r_max - r_th)/x_th * (x(101:end) - x_th);
x = x/x(end);
r_av = 6;
A = 1;
r = r_av + A * cos(2*pi*x);
% r(1:3*(n-1)/4+1) = r_av + A * cos(2*pi*x(1:3*(n-1)/4+1)/2.25);
% r(3*(n-1)/4+2:end) = r_av + A * cos(2*pi*(x(3*(n-1)/4+2:end)-0.75)/.75);
% taper_conv = pi/180 * 0.737;
% taper_div = pi/180 * 0.949;
% xp = 0.07 - 0.05 * sin(taper_conv);
% rp = 0.05265 - 0.05 * cos(taper_conv);
% xq = 0.07 + 0.05 * sin(taper_div);
% rq = 0.05265 - 0.05 * cos(taper_div);
% r_in = 0.0031;
% r_th = 0.00265;
% x_th = 0.07;
% r_fillet = 0.05;
% r_out = 0.0041;
% x_out = 0.245;
% for p = 1:length(x)
%     if x(p) < xp
%         r(p) = r_in + (rp - r_in)/xp * x(p);
%     elseif (x(p) >= xp) && (x(p) < xq)
%         r(p) = r_fillet + r_th - sqrt(r_fillet^2 - (x(p) - x_th)^2);
%     else
%         r(p) = rq + (r_out - rq)/(x_out - xq) * (x(p) - xq);
%     end
% end
r_exit = r(end);              % Exit diameter [m]
% x = x/x(end);
Length = x(end);
dx = Length/(n-1);

r_d = r/r(end) * r_exit; % Dimensional radius [m]
if circShape == true
    A_d = pi.*r_d.^2; % Area of the channel [m2]
    c_d = 2.*pi.*r_d; % Circumference of the channel [m]
else
    A_d = 2.*r_d.*crackLength; % Area of the channel [m2]
    c_d = 4.*r_d + 2.*crackLength; % Circumference of the channel [m]
end
% Find throat:
r_d_throat = min(r_d);
throat = find(r_d_throat == r_d); %n-loaction of throat
if length(throat) > 1
    disp("Warning: 2 locations for throat found")
end
throat = throat(1);
x_throat = x(throat);
% Make x, r and A non-dimensional
x_d = x .* L;
As = A_d(throat);
r = r_d./(As./L);
A = A_d/As;
c = c_d./(As./L); % Non-Dimensional circumference [-]
r_res = r(1);
r_throat = r(throat);
x_d_throat = x_d(throat);
% Generate initial flow profile with isentropic relations from 'minimoon.m'
[M_m, T_m, p_m, rho_m, V_m] = mini_moon(x, throat, A, Rg, kg);

% Dimensional
rho_inlet = rho_m(1);
p_inlet = p_m(1);
T_inlet = T_m(1);
M_inlet = M_m(1);
V_inlet = V_m(1);
rho_ref = rho_inlet;
p_ref = p_inlet;
T_ref = T_inlet;
M_ref = M_inlet;
% Make non-dimensional
T_m = T_m/T_inlet;
% if throat ~= n
%     nuc_period = throat + 0.01*n;
%     nuc_period = round(nuc_period);
%     T_m(throat:nuc_period) = T_m(throat) + (T_m(1)-T_m(throat))*(x(throat:nuc_period)-x(throat))./(x(nuc_period)-x(throat));
%     T_m(nuc_period:end) = T_m(1) - 0.1*(x(nuc_period:end)-x(nuc_period))./(x(end)-x(nuc_period));
% end
p_m = p_m/p_inlet;
rho_m = rho_m/rho_inlet;
exit = n;
tic
[~, V_d, p_d, ~, x_d, ~, T_d] = conservative_with_shock(x, dx, n, C, A, nt, kg, throat, Cx, Rg, pe, L, cv, T_ref, rho_ref, p_ref, r, c);
cons = toc; % end of run-time %
fprintf('Solution time for conservative form is %0.3g seconds', cons)

% plot([x_d_throat x_d_throat],[r_d(throat) -r_d(throat)],'k--')
% hold on
% plot(x_d,r_d,'k','linewidth',2)
% plot(x_d,-r_d,'k','linewidth',2)
% hold off
% D_exit = 2*r_d(end);
% D_throat = 2*r_d(throat);
% legend(['D_{exit}/D_{throat} = ', num2str(round(D_exit/D_throat,3))],['L = ', num2str(L), ' m'],['x_{throat}/L = ', num2str(round(x_throat,3))]);
% ylabel('Channel Radius [m]')
% xlabel('z [m]')

% T_cr = 647.096; % critical temperature of water [K]
% T_Wilson = min(T_d); % Wilson temperature
% Wilson_n = find(T_Wilson == T_d); % n-loaction of Wilson point
% Tsat_s0 = fzero(@(T) func(T, p_res, T_res, kg), T_res);
% Wilson = (Tsat_s0 - T_Wilson)/T_cr; % Wilson number
% t_act = simps(x_d(1:Wilson_n), 1./V_d(1:Wilson_n)); % activation time [s]
% Cr_av = Wilson/t_act; % average cooling rate [s^-1]

% Tsat_lg = linspace(273.16, 280, n);
% psat_lg = p_eq_lg(Tsat_lg);
% Tsat_sg = linspace(220, 273.16, n);
% psat_sg = p_eq_sg(Tsat_sg);
% A = [T_d; p_d];
% fid = fopen('15_pT.csv','w');
% fprintf(fid,'%12.8f,%12.8f\n',A);
% fclose(fid);
% fid = fopen('baseline_pT.csv','r');
% pT = fscanf(fid,'%g, %g',[2 Inf]);
% fclose(fid);
% T_d_base = pT(1,:);
% p_d_base = pT(2,:);
% fid = fopen('15_pT.csv','r');
% pT = fscanf(fid,'%g, %g',[2 Inf]);
% fclose(fid);
% T_d_15 = pT(1,:);
% p_d_15 = pT(2,:);
% hold on
% plot(T_d_base, p_d_base,'+')
% plot(T_d_15, p_d_15,'o')
% plot(T_d, p_d,'x')
% plot(273.16*ones(1,n), linspace(611.2,1e3,n),'k','linewidth',2)
% plot(Tsat_lg, psat_lg,'k','linewidth',2)
% plot(Tsat_sg, psat_sg,'k','linewidth',2)
% hold off
% xlabel('T [K]')
% ylabel('p [Pa]')
% legend({'Baseline','With wall friction and heat convection'},'Location','southeast')
% 
% function fun = func(T, p_res, T_res, kg)
%     fun = p_res * (T/T_res)^(kg/(kg-1)) - p_eq_sg(T);
% end
