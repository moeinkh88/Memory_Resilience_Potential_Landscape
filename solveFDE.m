function [LRdp,Lp,Hp,LamAttr,Q,Xall]=solveFDE(U,al,xt)

Nor=length(al);

TF=islocalmin(polyval(U,xt)); 
TF1=islocalmax(polyval(U,xt));
xmin=xt(TF); % critical points: bottom of valley
xmax=xt(TF1); % critical point: threshold 

X0=xmin(1)-(xmax-xmin(1)); % left state
MLX0=xmax-0.01; % middle state tend to left
MX0=xmax+0.01; % % middle state tend to right
RX0=xmin(2)+(xmin(2)-xmax); % rigth state

t0=0;
h=0.005;
F=@funPoly;
JF=@JfunPoly;
T=100;

true=0;
% epsilon=1e-4;
while true==0

% solve from 4 initial states to cover all states for potential landscape
for J=1:Nor
%          [t,x] = FDE_PI12_PC(al(J),F,t0,T,X0,h); 
                  [t,x] = FDE_PI2_IM(al(J),F,JF,t0,T,X0,h); 
    X(J,:)=x;
%          [MLt,MLx] = FDE_PI12_PC(al(J),F,t0,T,MLX0,h); 
                  [MLt,MLx] = FDE_PI2_IM(al(J),F,JF,t0,T,MLX0,h); 
    MLX(J,:)=MLx;
%          [Mt,Mx] = FDE_PI12_PC(al(J),F,t0,T,MX0,h); 
                  [Mt,Mx] = FDE_PI2_IM(al(J),F,JF,t0,T,MX0,h); 
    MX(J,:)=Mx;
%          [Rt,Rx] = FDE_PI12_PC(al(J),F,t0,T,RX0,h); 
                  [Rt,Rx] = FDE_PI2_IM(al(J),F,JF,t0,T,RX0,h); 
    RX(J,:)=Rx;
end

if abs(x(end)-xmin(1))<2*h && ...
        abs(MLx(end)-xmin(1))<2*h && ...
        abs(Mx(end)-xmin(2))<2*h && ...
        abs(Rx(end)-xmin(2))<2*h
    
    true=1;
else
    T=T+100;
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
NX=.5*(RX0-X0)/h; %number of points in the interpolation interval

for J=1:Nor
[C,ia] = unique(Xall1(:,J)); % sort unique
Xall(:,J)=linspace(X0,RX0, NX); % interval of interpolation
DX(:,J) = interp1(C,DX1(ia,J),Xall(:,J),'pchip'); % interpilation

Q(:,J) = -cumtrapz(Xall(:,J),DX(:,J)); % potential values

TF = islocalmin(Q(:,J)); %find attraction states
MinU(:,J)=TF;  
 
TF1 = islocalmax(Q(:,J));
MaxU(:,J)=TF1;

Hp(J)=Q(MaxU(:,J),J); % height of potantial landscape on threshold (peak)
Lp(J,:)=Q(MinU(:,J),J);  % height of potantial landscape in attractions
LRdp(J,:)=abs(Hp(J)-Lp(J,:)); % dept of basin of attractions (left and right)

H=((C(end)-C(1))/NX);
DQ(:,J)=diff(Q(:,J))/H;
indx=Xall(TF,J);
    XLattr=find(Xall(:,J)==indx(1));
    XRattr=find(Xall(:,J)==indx(2));
LamAttr(J,1)=-(DQ(XLattr+1,J)-DQ(XLattr,J))/H;
LamAttr(J,2)=-(DQ(XRattr+1,J)-DQ(XRattr,J))/H;
end

end