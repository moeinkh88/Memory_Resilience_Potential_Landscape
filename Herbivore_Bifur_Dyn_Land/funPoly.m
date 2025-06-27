function dx=funPoly(t,x)

global par

a1=par(3);
a2=par(2);
a3=par(1);

dx= a3*x^3+a2*x^2+a1*x;
end