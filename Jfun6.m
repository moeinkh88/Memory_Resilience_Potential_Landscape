function dx=Jfun6(t,x)
global rr KK lambda aa 

dx= rr*(1-x/KK)-rr/KK*x-lambda/(x+aa)+lambda*x/(x+aa)^2;

end