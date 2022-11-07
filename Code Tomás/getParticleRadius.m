function R = getParticleRadius(dR_dx, x_d)
% Set particle radius to zero @z=0       
R(1) = 0;
% Update maximal particle radius, Change to Simpson's rule %
for i = 2:length(dR_dx)
%    R(i) = dR_dx(i-1)*dx_d + R(i-1);
    R(i) = simps(x_d(1:i),dR_dx(1:i));
%     R(i) = trapz(x_d(1:i),dR_dx(1:i));
    if R(i) < 0 % If particle radius < 0 , set as zero 
        R(i) = 0;
        dR_dx(i) = 0;
        if i == 2
            dR_dx(1) = 0;
        end
                
    end
end
end