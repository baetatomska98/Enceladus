function [M, Tstat, pstat, rho, V] = mini_moon(x, throat, A, R, k)  
Tres = 273.16;
Pres = 611.2;
% rhores = Pres/(R*Tres);

M = zeros(1, length(x));
Tr = zeros(1, length(x));
pr = zeros(1, length(x));
% rhor = zeros(1, length(x));
for i = 1:length(x)
    if i < throat
        [M(i), Tr(i), pr(i), ~, ~] = flowisentropic(k, A(i), 'sub');
    else
        [M(i), Tr(i), pr(i), ~, ~] = flowisentropic(k, A(i), 'sup');
    end
end
Tstat = Tr * Tres;
pstat = pr * Pres;
rho = pstat./(R*Tstat); %rhor * rhores;
V = M .* sqrt(k*R*Tstat);
end