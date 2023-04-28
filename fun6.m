function dx=fun6(t,x)
global rr KK lambda aa 

dx= (rr*x*(1-x/KK)-lambda*x/(x+aa));
end