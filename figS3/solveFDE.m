function [Q,Xall]=solveFDE(par,al)

xt=0:.00001:2.5; % steps for polynomial function

Nor=length(al);

X0=0; % left state

RX0=xt(end);

t0=0;
h=0.005;
F=@funPoly;
JF=@JfunPoly;
T=1000;

true=0;
% epsilon=1e-4;
while true==0

% solve from 4 initial states to cover all states for potential landscape
for J=1:Nor

                  [Rt,Rx] = FDE_PI2_IM(al(J),F,JF,t0,T,RX0,h,par); 
    RX(J,:)=Rx;
end

if        abs(Rx(end)-X0)<2*h
    
    true=1;
else
    T=T+1000;
    clear X MLX MX RX
end
end

% Rate of change; right side of the main equation
dxR=diff(RX')./h;

J=3:length(RX)-2;
dxR(J-1,:)=   1/(12*h).*(RX(:,J-2)' - 8.*RX(:,J-1)' + 8.*RX(:,J+1)' - RX(:,J+2)');

% Sorting the states in the increasing order
DX1=flip(dxR);
Xall1=flip(RX(:,1:end-1)');

% interpolation the generated data
NX=20*(RX0-X0)/h; %number of points in the interpolation interval

for J=1:Nor
% [C,ia] = unique(Xall1(:,J)); % sort unique
% Xall(:,J)=linspace(X0,RX0, NX); % interval of interpolation
% DX(:,J) = interp1(C,DX1(ia,J),Xall(:,J),'pchip'); % interpilation
Q(:,J) = -cumtrapz(Xall1(:,J),DX1(:,J)); % potential values
Xall=Xall1;
end

end