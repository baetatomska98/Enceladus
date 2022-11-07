% cdiff - calculates central difference of a input function

function f1 = cdiff(f,i)

f1 = (f(i+1) - 2*f(i) + f(i-1));

end