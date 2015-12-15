%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Newton's Method           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = 1.4;
M = 0:.01:4;
A = 1./M.*((2./(gamma+1)).*(1+(gamma-1)./2.*M.^2)).^((gamma+1)./...
    (2.*(gamma-1)));
%  This is the algebraic equation for A - a function of M

Newton(

syms m
function f = func1()
f = sym(1./M.*((2./(gamma+1)).*(1+(gamma-1)./2.*M.^2)).^((gamma+1)./...
    (2.*(gamma-1)))-A);
end 
f2 = subs(fun1,m);
%Creates symbolic function that is used to take the derivative of A

d = diff(f2,m);
f3 = vectorize(d);
n = 1;
x(1) = 0;
error = 1;

while error > .0001
    x(n+1) = x(n)-f(x(n))/f3(x(n));
    error = abs(x(n)-x(n+1));
    n = n+1;
end


hFig = figure(1);
plot(M,A)

%  This outputs a plot of the corresponding function
    