function [Xall,Q]=PotenDepth(par,al)

% X0=0; % left state
% MLX0=.8432-0.001; % middle state tend to left
% MX0=.8432+0.001; % % middle state tend to right
% RX0=3; % right state
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
% RX0=xt(Indx2(1)+indx); % right state
RX0=2.5;

t0=0;
h=0.001;
F=@funPoly;
JF=@JfunPoly;
T=100;

true=0;
% epsilon=1e-4;
while true==0

% solve from 4 initial states to cover all states for potential landscape
for J=1:1
%          [t,x] = FDE_PI12_PC(al(J),F,t0,T,X0,h); 
                  [t,x] = FDE_PI2_IM(al(J),F,JF,t0,1,X0,h,par); 
    X(J,:)=x;
%          [MLt,MLx] = FDE_PI12_PC(al(J),F,t0,T,MLX0,h); 
                  [MLt,MLx] = FDE_PI2_IM(al(J),F,JF,t0,T,MLX0,h,par); 
    MLX(J,:)=MLx;
%          [Mt,Mx] = FDE_PI12_PC(al(J),F,t0,T,MX0,h); 
                  [Mt,Mx] = FDE_PI2_IM(al(J),F,JF,t0,T,MX0,h,par); 
    MX(J,:)=Mx;
%          [Rt,Rx] = FDE_PI12_PC(al(J),F,t0,T,RX0,h); 
                  [Rt,Rx] = FDE_PI2_IM(al(J),F,JF,t0,T,RX0,h,par); 
    RX(J,:)=Rx;
end

if      abs(MLx(end)-0)<2*h && ...
        abs(Mx(end)-xmin)<2*h && ...
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
DX(:,J) = interp1(C,DX1(ia,J),Xall(:,J),'pchip'); % interpolation

Q(:,J) = -cumtrapz(Xall(:,J),DX(:,J)); % potential values

% TF = islocalmin(Q(:,J)); %find attraction states
% 
%  
% TF1 = islocalmax(Q(:,J));
% MaxU(:,J)=TF1;
% 
% Hp(J)=Q(TF1,J); % height of potential landscape on threshold (peak)
% Lp(J,1)=0;
% Lp(J,2)=Q(TF,J);  % height of potential landscape in attractions
% LRdp(J,:)=abs(Hp(J)-Lp(J,:)); % dept of basin of attractions (left and right)
% 
% H=((C(end)-C(1))/NX);
% DQ(:,J)=diff(Q(:,J))/H;
% indx=Xall(TF,J);
% %     XLattr=find(Xall(:,J)==indx(1));
%     XRattr=find(Xall(:,J)==indx);
% LamAttr(J,1)=-(DQ(2,J)-DQ(1,J))/H;
% LamAttr(J,2)=-(DQ(XRattr+1,J)-DQ(XRattr,J))/H;
end

end