function dx=fun5(t,x)
global rr KK lambda aa 

dx= (rr*x*(1-x/KK)-lambda*x^2/(x^2+aa^2));
end