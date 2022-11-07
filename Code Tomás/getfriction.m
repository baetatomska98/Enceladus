% friction calc.
function [friction] = getfriction(tau, c, dx)
%friction = tau.*pi.*2.*r;%.*L./sqrt(A_d(throat)); % the "L./sqrt(A_d(throat)" is to non-dimensionalize 
friction = tau.*c.*dx;
end