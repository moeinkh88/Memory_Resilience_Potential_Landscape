function dx=funC(t,x,P)
global r K A B

p=P(1);pt0=P(2);
%%%pertubation
if t>pt0
b=p+B;
else
    b=B;
end
%%%

dx= (r*x*(1-x/K)-b*x/(x+A));
end