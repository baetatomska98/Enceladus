function [fp, gn] = solid_frac_MC(Q,dx,T,rho,R_tab,dR_dx_tab,A_tab,Rg)
rho_grain = 920; %density of ice [kg/m3]
%[gn_tab S] = gamma_nuc_MC(T,rho); %[/cm3s] Only for verification of C100 & C101
gn_tab = getNucleationRate(T, rho, Rg); %[/cm3s]
gn_tab = gn_tab.*10^6; % [/m^3s] 
for i = 1:length(A_tab) % Computes whole length %
    A = A_tab(1:i);
    dR_dx = dR_dx_tab(1:i);
    R = R_tab(1:i);
    gn = gn_tab(1:i);
    G = gn.*(R(end)-R).^2.*dR_dx(end).*A.*heaviside(R(end)-R);
    if i == 1
        continue
    end
    if length(dx) > 1
        fp(i) = (4*pi*rho_grain/Q)*dx(i-1)*simps(G);
        %fp(i) = (4*pi*rho_grain/Q)*dx(i-1)*trapz(G);
    else
        fp(i) = (4*pi*rho_grain/Q)*dx*simps(G);
        %fp(i) = (4*pi*rho_grain/Q)*dx*trapz(G);
    end
    
end 
end