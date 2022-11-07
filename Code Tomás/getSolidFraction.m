function f = getSolidFraction(fp, x_d)
% Set solid fraction and its derivative to zero at z=0
f(1) = 0;
fp(1) = 0;
for i = 2:length(fp)
    %f(i) = f(i-1) + fp(i-1)*dx_d;
    f(i) = simps(x_d(1:i),fp(1:i));
    %f(i) = trapz(x_d(1:i),fp(1:i));
    if f(i)>=1 % If solid-fraction > 1, set as 1  
        f(i) = 1;
        fp(i) = 0;
    elseif f(i) <= 0 
        f(i) = 0; % If solid-fraction <= 0, set as 0 
    end
end
end