function dx=Jfun(t,x,b)
global r K A 

dx= r*(1-x/K)-r/K*x-b/(x+A)+b*x/(x+A)^2;

end