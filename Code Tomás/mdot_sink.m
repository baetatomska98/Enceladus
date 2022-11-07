% Deposited Mass %
function [m_deposition] = mdot_sink(f, rho, A, V)
%m_deposition = (f(i)-f(i-1)) * rho(i) * A(i) * V(i); % gas -> solid "deposition" SINK TERM
%m_deposition = fp(i)* dx_d * rho(i) * A(i) * V(i);
%m_deposition = fp(i)* dx_d * A(i);% * V(i);
% add additional terms if needed
%m_deposition = fp .* dx_d .* rho .* A .* V;
m_deposition = f .* rho .* A .* V; % Non-dimensional
end