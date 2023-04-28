function dx=Jfun(t,x)
global rr KK lambda aa 

dx= (rr*(1-x/KK)-rr/KK*x -2*lambda*x/(x^2+aa^2)+2*lambda*x^3/((x^2+aa^2)^2));

end