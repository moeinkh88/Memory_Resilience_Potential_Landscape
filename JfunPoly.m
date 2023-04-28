function dx=JfunPoly(t,x)
global par

a1=par(3);
a2=par(2);
a3=par(1);

dx= 3*a3*x^2+2*a2*x+a1;
end