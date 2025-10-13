function [Xall,Q]=PotenDepth(par,A,al)

% X0=0; % left state
% MLX0=.8432-0.001; % middle state tend to left
% MX0=.8432+0.001; % % middle state tend to right
% RX0=3; % rigth state
% xmin=1.9568;

xt=0:.00001:2.5; % steps for polynomial function
U=-polyint(par);
U1=polyval(U,xt);
TF=islocalmin(U1);
TF1=islocalmax(U1);

Nor=length(al);

xmin=xt(TF); % critical points: bottom of valley
xmax=xt(TF1); % critical point: threshold 

indx=find(TF==1);
ind=0; epsi1=0.0001;epsi2=0.00001;
while ind==0

Indx2=find((abs(U1(indx:end)-U1(TF1))<epsi2));

if isempty(Indx2)==1
    epsi2=epsi2*2;
else
    ind=1;
end
end

X0=0; % left state
MLX0=xmax-0.001; % middle state tend to left
MX0=xmax+0.001; % % middle state tend to right
% RX0=xt(Indx2(1)+indx); % rigth state
RX0=2.5;

t0=0;
h=0.001;
F=@funPoly_den; % Herbivory function converted to polynomial
JF=@JfunPoly_den; % Jacobian of polynomial function
T=100;

true=0;
% epsilon=1e-4;
while true==0

% solve from 4 initial states to cover all states for potential landscape
for J=1:1
%          [t,x] = FDE_PI12_PC(al(J),F,t0,T,X0,h); 
                  [t,x] = FDE_PI2_IM(al(J),F,JF,t0,1,X0,h,[par,A]); 
    X(J,:)=x;
%          [MLt,MLx] = FDE_PI12_PC(al(J),F,t0,T,MLX0,h); 
                  [MLt,MLx] = FDE_PI2_IM(al(J),F,JF,t0,T,MLX0,h,[par,A]); 
    MLX(J,:)=MLx;
%          [Mt,Mx] = FDE_PI12_PC(al(J),F,t0,T,MX0,h); 
                  [Mt,Mx] = FDE_PI2_IM(al(J),F,JF,t0,T,MX0,h,[par,A]); 
    MX(J,:)=Mx;
%          [Rt,Rx] = FDE_PI12_PC(al(J),F,t0,T,RX0,h); 
                  [Rt,Rx] = FDE_PI2_IM(al(J),F,JF,t0,T,RX0,h,[par,A]); 
    RX(J,:)=Rx;
end

if      abs(MLx(end)-0)<2*h && ...
        abs(Mx(end)-xmin)<4*h && ...
        abs(Rx(end)-xmin)<2*h
    
    true=1;
else
    T=T+200;
    clear X MLX MX RX
end
end

% Rate of change; right side of the main equation
dx=diff(X')./h;
dxML=diff(MLX')./h;
dxM=diff(MX')./h;
dxR=diff(RX')./h;

J=3:length(X)-2;
dx(J-1,:)=   1/(12*h).*(X(:,J-2)' - 8.*X(:,J-1)' + 8.*X(:,J+1)' - X(:,J+2)');
dxML(J-1,:)=   1/(12*h).*(MLX(:,J-2)' - 8.*MLX(:,J-1)' + 8.*MLX(:,J+1)' - MLX(:,J+2)');
dxM(J-1,:)=   1/(12*h).*(MX(:,J-2)' - 8.*MX(:,J-1)' + 8.*MX(:,J+1)' - MX(:,J+2)');
dxR(J-1,:)=   1/(12*h).*(RX(:,J-2)' - 8.*RX(:,J-1)' + 8.*RX(:,J+1)' - RX(:,J+2)');

% Sorting the states in the increasing order
DX1=cat(1,dx,flip(dxML),dxM,flip(dxR));
Xall1=cat(1,X(:,1:end-1)',flip(MLX(:,1:end-1)'),MX(:,1:end-1)',flip(RX(:,1:end-1)'));

% interpolation the generated data
NX=2*(RX0-X0)/h; %number of points in the interpolation interval

for J=1:1
[C,ia] = unique(Xall1(:,J)); % sort unique
Xall(:,J)=linspace(X0,RX0, NX); % interval of interpolation
DX(:,J) = interp1(C,DX1(ia,J),Xall(:,J),'pchip'); % interpilation

Q(:,J) = -cumtrapz(Xall(:,J),DX(:,J)); % potential values

end

%%
function dx=funPoly_den(t,x,par)


a1=par(3);
a2=par(2);
a3=par(1);
A=par(4);

dx= (a3*x^3+a2*x^2+a1*x)/(x+A);
end


%%
function dx=JfunPoly_den(t,x,par)

a1=par(3);
a2=par(2);
a3=par(1);
A=par(4);

num = 2*a3*x.^3 + (3*a3*A + a2)*x.^2 + 2*a2*A*x + a1*A;
    den = (A + x).^2;
    dx   = num ./ den;
end

end