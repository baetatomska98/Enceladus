% Latent heat as a source term for the energy eq. [J/s]%
function Edot_L = Edot_source(rho, f, fp, dx_d, V, L_h, A, cv, Rg, dx, fp_nd, Lh_nd, x, T, y)
%Edot_L = rho .* f .* A .* V .* (Lh_nd + (y - 1)*cv/Rg .* T);
Edot_L = rho .* f .* A .* V .* Lh_nd;
%Edot_L = rho(1) .* f .* A .* V .* Lh_nd; % Multiply with rho(1) = 1
%Edot_L = rho(1) .* fp_nd .* dx .* A .* V .* Lh_nd; % Multiply with rho(1) = 1
%Edot_L = rho .* fp_nd .* dx .* A .* V .* Lh_nd;
end