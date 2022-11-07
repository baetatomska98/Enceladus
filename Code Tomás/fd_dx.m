% fd_dx - calculates forward difference of a input function

function f1 = fd_dx(f,dx,i)
if length(dx)>1
    f1 = (f(i+1) - f(i))/dx(i);
else
    f1 = (f(i+1) - f(i))/dx;
end
end