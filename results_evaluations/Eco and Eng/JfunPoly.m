function dx=JfunPoly(t,x,par)

pt0=par(5);
pt1=par(6);
p=par(4);
a1=par(3);
a2=par(2);
a3=par(1);


%%%pertubation
if t>pt0 && t<pt1
a1=par(3)-p;
end
%%%

dx= 3*a3*x^2+2*a2*x+a1;
end