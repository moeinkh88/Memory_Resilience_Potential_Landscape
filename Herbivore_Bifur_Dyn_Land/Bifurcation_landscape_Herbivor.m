clc
clear 

global r K A B

% coefficients
r=.8;
K=6;
A=.1;
B=1;

b=.4:.01:2;

%%
syms c xx

eqn=(r*xx*(1-xx/K)-c*xx/(xx+A)); % equation function herbivore model
sc=solve(eqn,xx) % to find the equilibrium points

% c1=0.16:.001:1;
% c2=0.16:.001:.682;
% c3=0:.001:.682;
c1=0.08:.001:2;
c2=0.08:.001:1.24 ;
c3=0:.001:1.24;
Eq1=                          0.*c1;
Eq2= 59/20 - (3721 - 3000*c2).^(1/2)/20;
Eq3= (3721 - 3000*c3).^(1/2)/20 + 59/20;
figure
p=plot(c1,Eq1,'k',c2,Eq2,'k--',c3,Eq3,'k'); % Bifurcation diagram

set(p,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Equilibrium density')
xlabel('Paramter B')

%% Converting equation to a polynomial format 
h=0.01;
x=0:h:6;

a1=A*r-B;
a2=r*(1-A/K);
a3=-r/K;

%% Evaluating potential landscape

global par
xt=0:.00001:9; % steps for polynomial function
al=1:-.2:.8; % order of derivatives (1-memory) 
Nor=length(al);

par=[a3 a2 a1];

U=-polyint(par); % potential of polynomial fun; -int(poly)
U1=polyval(U,xt);
TF=islocalmin(U1);
TF1=islocalmax(U1);

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
MLX0=xmax-0.001; % unstable state tending to left
MX0=xmax+0.001; % % unstable state tending to right
RX0=xt(Indx2(1)+indx); % right state

t0=0;
h=0.005;
F=@funPoly; % Herbivory function converted to polynomial
JF=@JfunPoly; % Jacobian of polynomial function
T=300; % final time for the dynamics

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

if         abs(MLx(end))<2*h && ...
        abs(Mx(end)-xmin(1))<2*h && ...
        abs(Rx(end)-xmin(1))<2*h
    
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

J=3:length(MX)-2;
dx(J-1,:)=   1/(12*h).*(X(:,J-2)' - 8.*X(:,J-1)' + 8.*X(:,J+1)' - X(:,J+2)');
dxML(J-1,:)=   1/(12*h).*(MLX(:,J-2)' - 8.*MLX(:,J-1)' + 8.*MLX(:,J+1)' - MLX(:,J+2)');
dxM(J-1,:)=   1/(12*h).*(MX(:,J-2)' - 8.*MX(:,J-1)' + 8.*MX(:,J+1)' - MX(:,J+2)');
dxR(J-1,:)=   1/(12*h).*(RX(:,J-2)' - 8.*RX(:,J-1)' + 8.*RX(:,J+1)' - RX(:,J+2)');

% Sorting the states in the increasing order
DX1=cat(1,dx,flip(dxML),dxM,flip(dxR));
Xall1=cat(1,X(:,1:end-1)',flip(MLX(:,1:end-1)'),MX(:,1:end-1)',flip(RX(:,1:end-1)'));

% interpolation the generated data
NX=2*(RX0-0)/h; %number of points in the interpolation interval

for J=1:Nor
[C,ia] = unique(Xall1(:,J)); % sort unique
Xall(:,J)=linspace(0,RX0, NX); % interval of interpolation
DX(:,J) = interp1(C,DX1(ia,J),Xall(:,J),'pchip'); % interpolation

Q(:,J) = -cumtrapz(Xall(:,J),DX(:,J)); % potential values

TF = islocalmin(Q(:,J)); %find attraction states
MinU(:,J)=TF;  
 
TF1 = islocalmax(Q(:,J));
MaxU(:,J)=TF1;

Hp(J)=Q(MaxU(:,J),J); % height of potential landscape on threshold (peak)
Lp(J,:)=Q(MinU(:,J),J);  % height of potential landscape in attractions
LRdp(J,:)=abs(Hp(J)-Lp(J,:)); % dept of basin of attractions (left and right)

H=((C(end)-C(1))/NX);
DQ(:,J)=diff(Q(:,J))/H;
indx=Xall(TF,J);
    XRattr=find(Xall(:,J)==indx);
LamAttr(J)=-(DQ(XRattr+1,J)-DQ(XRattr,J))/H;
end


%% plotting 
% colors=[0.4940 0.1840 0.5560
% 0.4660 0.6740 0.1880
% 0.6350 0.0780 0.1840];
colors=[0 0.4470 0.7410
    0.8500 0.3250 0.0980];

figure 
subplot(2,2,1) % bifurcation diagram of equilibrium density
p=plot(c1,Eq1,'k',c2,Eq2,'k--',c3,Eq3,'k');

set(p,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Equilibrium density')
xlabel('Paramter B')

hold on
xline(1,':')
text(1,xmin+.2,'X_{s2}','FontSize',14)
text(1-.2,xmax,'X_{u}','FontSize',14)
text(1,0+.4,'X_{s1}','FontSize',14)

subplot(2,2,2) % System dynamics from 3 initial conditions; 2 close the unstable points MX, MLX, and one from the right side of the right valley, RX

for i=1:length(al)

hold on

% p1=plot(t,X(i,:),'color',colors(i,:));
p1ML=plot(MLt,MLX(i,:),'color',colors(i,:));

p1M=plot(Mt,MX(i,:),'color',colors(i,:));
% legend('Memoryless','Memory=0.2')
% set(get(get(p1R,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

p1R=plot(Rt,RX(i,:),'color',colors(i,:));
legend('Memory=0','Memory=0.2')
set(get(get(p1R,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p1M,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(p1ML,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


% set(p1,'LineWidth',3)
% set(gca,'FontSize',14)
set(p1ML,'LineWidth',3)
set(gca,'FontSize',14)
set(p1M,'LineWidth',3)
set(gca,'FontSize',14)
set(p1R,'LineWidth',3)
set(gca,'FontSize',14)
end

xlabel('time')
ylabel('States')
axis([0,20,0,RX0])

text(1,xmin-.4,'X_{s2}','FontSize',14)
text(15,xmax,'X_{u}','FontSize',14)
text(1,0+.4,'X_{s1}','FontSize',14)
hold on
p1R=yline(xmax,':');
p1M=yline(xmin,':');
set(get(get(p1R,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p1M,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

subplot(2,2,3) % Potential landscape 
hold on
for i=1:2
   p4=plot(Xall(:,i),Q(:,i),'color',colors(i,:));
hold on

set(p4,'LineWidth',3)
set(gca,'FontSize',14) 
end
xlabel('States')
ylabel('Potential energy')
axis([0, 5.5, Q(MinU(:,1),1),Q(MaxU(:,1),1)])

subplot(2,2,4) % slope of potential landscape 
hold on
for i=1:2
   p4=plot(Xall(:,i),DX(:,i),'color',colors(i,:));
hold on

set(p4,'LineWidth',3)
set(gca,'FontSize',14) 
end
xlabel('States')
ylabel('Growth rates')
yline(0,':')
axis([0 5 -1.5 .8 ])
