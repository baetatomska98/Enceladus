function [fp_nd, Lh_nd, Rg_nd, E_nd] = nondimensionalize(fp, L, L_h, cv, T_ref, Rg, E, rho_ref, kg)
fp_nd = fp .* L;
Lh_nd = L_h/(Rg*T_ref);
Rg_nd = Rg/cv;
a0 = sqrt(Rg*T_ref*kg);
%E_nd = E.*(sqrt(y)./(rho_ref.*a0)); % Approach 1
E_nd = E./(rho_ref.*a0); % Approach 2
end