function dx=funC(t,x)
global rr KK aa C

dx= (rr*x*(1-x/KK)-C*x/(x+aa));
end