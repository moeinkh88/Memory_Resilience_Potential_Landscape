function dx=fun(t,x,b)
global r K A  

dx= (r*x*(1-x/K)-b*x/(x+A));
end