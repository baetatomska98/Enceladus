function p = p_eq_lg(T)
% Water liquid - gas equilibrium pressure
p = 610.8.*exp(-5.1421.*log(T./273.15)-6828.77*(1./T-1/273.15));