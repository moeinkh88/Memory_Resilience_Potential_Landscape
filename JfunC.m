function dx=JfunC(t,x)
global rr KK aa C

dx= rr*(1-x/KK)-rr/KK*x-C/(x+aa)+C*x/(x+aa)^2;

end