function [rho_d, T_d, p_d, V_d, A_d, x_d, dx_d, r_d] = dimensionalize(rho, T, p, M, kg, Rg, x, dx, rho_ref, T_ref, p_ref, L, r, As, circShape, crackLength)
rho_d = rho * rho_ref;
T_d = T * T_ref;
p_d = p * p_ref;
V_d = M .* sqrt(kg * Rg .* T_d);
r_d = r.*(As./L);
if circShape == false
    A_d = 2.*r_d.*crackLength;
else
    A_d = pi.*r_d.^2;
end
x_d = x .* L;
dx_d = dx .* L;
end