% rd_dx - calculates rearward differences of the input function

function f1 = rd_dx(f,dx,i)
if length(dx)>1
    f1 = (f(i) - f(i-1))/dx(i-1); %dx(i-1) ?
else
    f1 = (f(i) - f(i-1))/dx;
end
end