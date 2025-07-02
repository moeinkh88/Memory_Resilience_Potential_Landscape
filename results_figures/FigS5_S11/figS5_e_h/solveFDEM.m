function [QR,XallR]=solveFDEM(par,al)

xt=0:.00001:1.9568; % steps for polynomial function

Nor=length(al);

X0=xt(1); % left state

RX0=xt(end);

t0=0;
h=0.005;
F=@funCP;
JF=@JfunCP;
T=1000;

true=0;
% epsilon=1e-4;
while true==0

% solve from 4 initial states to cover all states for potential landscape
for J=1:Nor

                  [Rt,Rx] = FDE_PI2_IM(al(J),F,JF,t0,T,RX0,h,par); 
    RX(J,:)=Rx;
end

if par(1)>.275 && abs(Rx(end)-X0)<2*h        
        true=1;
elseif par(1)<.2755 && abs(Rx(end)-RX0)<2*h        
    true=1;
else
        T=T+1000;
        clear X MLX MX RX
end
end

% Rate of change; right side of the main equation
%before perturb
RXL=RX(1,1:20*200);
RXR=RX(1,20*200:end);

dxRL=diff(RXL')./h;
J=3:length(RXL)-2;
dxRL(J-1,:)=   1/(12*h).*(RXL(:,J-2)' - 8.*RXL(:,J-1)' + 8.*RXL(:,J+1)' - RXL(:,J+2)');

dxRR=diff(RXR')./h;
J=3:length(RXR)-2;
dxRR(J-1,:)=   1/(12*h).*(RXR(:,J-2)' - 8.*RXR(:,J-1)' + 8.*RXR(:,J+1)' - RXR(:,J+2)');

% Sorting the states in the increasing order
DX1L=flip(dxRL);
Xall1L=flip(RXL(:,1:end-1)');

DX1R=flip(dxRR);
Xall1R=flip(RXR(:,1:end-1)');

% interpolation the generated data
NX=20*(RX0-X0)/h; %number of points in the interpolation interval

for J=1:Nor
% [C,ia] = unique(Xall1(:,J)); % sort unique
% Xall(:,J)=linspace(X0,RX0, NX); % interval of interpolation
% DX(:,J) = interp1(C,DX1(ia,J),Xall(:,J),'pchip'); % interpilation
QL(:,J) = -cumtrapz(Xall1L(:,J),DX1L(:,J)); % potential values
XallL=Xall1L;

QR(:,J) = -cumtrapz(Xall1R(:,J),DX1R(:,J)); % potential values
XallR=Xall1R;
end

end