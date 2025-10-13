clc
clear 

%% Orders of the fractional derivative
al = 1:-.2:.8;     % → [1  0.8]
Nor = numel(al);

xt = 0:0.00001:9;       % still outside the loop – doesn’t depend on params
F=@funPoly; % Herbivory function converted to polynomial
JF=@JfunPoly; % Jacobian of polynomial function
%% ------------------------------------------------------------------
%  Parameter sets
%  ------------------------------------------------------------------
%  α = 1   (no‑memory)   → new values you gave
r_nomem = 0.739859;
K_nomem = 5.88714;
A_nomem = 0.0688134;
B_nomem = 0.944932;

%  α ≠ 1   (memory)      → keep your old values
r_mem = 0.8;
K_mem = 6;
A_mem = 0.1;
B_mem = 1;


b=.4:.01:2;


h=0.01;
x=0:h:6;


MLX1 = cell(Nor,1);
MX1  = cell(Nor,1);
RX1  = cell(Nor,1);
MLt1 = cell(Nor,1);
Mt1 = cell(Nor,1);
Rt1 = cell(Nor,1);
Q1 = cell(Nor,1);
Dx1 = cell(Nor,1);
Xall2 = cell(Nor,1);

%% Evaluating potential landscape
for j = 1:Nor
    
    % ---------- 1. select the right parameter set -------------------
    if abs(al(j) - 1) < eps          % α = 1 → no memory
        r = r_nomem;   K = K_nomem;
        A = A_nomem;   B = B_nomem;
    else                            % α ≠ 1 → memory
        r = r_mem;     K = K_mem;
        A = A_mem;     B = B_mem;
    end
    
    % ---------- 2. every quantity that depends on r,K,A,B ------------
    a1 = A*r - B;
    a2 = r*(1 - A/K);
    a3 = -r/K;
    par = [a3 a2 a1];

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
MLX0=xmax-0.002; % unstable state tending to left
MX0=xmax+0.002; % % unstable state tending to right
RX0=xt(Indx2(1)+indx); % right state

t0=0;
h=0.005;
T=300; % final time for the dynamics

true=0;
% epsilon=1e-4;
while true==0

% solve from 4 initial states to cover all states for potential landscape

                  [t,X] = FDE_PI2_IM(al(j),F,JF,t0,T,X0,h,par); 
    
                  [MLt,MLX] = FDE_PI2_IM(al(j),F,JF,t0,T,MLX0,h,par); 
    
                  [Mt,MX] = FDE_PI2_IM(al(j),F,JF,t0,T,MX0,h,par); 
    
                  [Rt,RX] = FDE_PI2_IM(al(j),F,JF,t0,T,RX0,h,par); 
    


if         abs(MLX(end))<2*h && ...
        abs(MX(end)-xmin(1))<2*h && ...
        abs(RX(end)-xmin(1))<2*h

    
    true=1;
else
    T=T+100;
    clear X(j,:) MLX(j,:) MX(j,:) RX(j,:)
end
end

% Rate of change; right side of the main equation
dx=diff(X')./h;
dxML=diff(MLX')./h;
dxM=diff(MX')./h;
dxR=diff(RX')./h;

J=3:length(MX)-2;
dx(J-1)=   1/(12*h).*(X(J-2)' - 8.*X(J-1)' + 8.*X(J+1)' - X(J+2)');
dxML(J-1)=   1/(12*h).*(MLX(J-2)' - 8.*MLX(J-1)' + 8.*MLX(J+1)' - MLX(J+2)');
dxM(J-1)=   1/(12*h).*(MX(J-2)' - 8.*MX(J-1)' + 8.*MX(J+1)' - MX(J+2)');
dxR(J-1)=   1/(12*h).*(RX(J-2)' - 8.*RX(J-1)' + 8.*RX(J+1)' - RX(J+2)');

% Sorting the states in the increasing order
DX1=cat(1,dx,flip(dxML),dxM,flip(dxR));
Xall1=cat(1,X(1:end-1)',flip(MLX(1:end-1)'),MX(1:end-1)',flip(RX(1:end-1)'));

% interpolation the generated data
NX=2*(RX0-0)/h; %number of points in the interpolation interval

[C,ia] = unique(Xall1); % sort unique
Xall=linspace(0,RX0, NX); % interval of interpolation
DX = interp1(C,DX1(ia),Xall,'pchip'); % interpolation

Q= -cumtrapz(Xall,DX); % potential values

TF = islocalmin(Q); %find attraction states
MinU=TF;  
 
TF1 = islocalmax(Q);
MaxU=TF1;

H=((C(end)-C(1))/NX);
DQ=diff(Q)/H;
indx=Xall(TF);
    XRattr=find(Xall==indx);
LamAttr(J)=-(DQ(XRattr+1)-DQ(XRattr))/H;

MLX1{j}=MLX;
MLt1{j}=MLt;
MX1{j}=MX;
Mt1{j}=Mt;
RX1{j}=RX;
Rt1{j}=Rt;
Q1{j}=Q;
Dx1{j}=DX;
Xall2{j}=Xall;
end


%% plotting 
% colors=[0.4940 0.1840 0.5560
% 0.4660 0.6740 0.1880
% 0.6350 0.0780 0.1840];
colors=[0 0.4470 0.7410
    0.8500 0.3250 0.0980];

figure 
subplot(2,2,1) % bifurcation diagram of equilibrium density

subplot(2,2,2) % System dynamics from 3 initial conditions; 2 close the unstable points MX, MLX, and one from the right side of the right valley, RX

for i=1:length(al)

hold on

% p1=plot(t,X(i,:),'color',colors(i,:));
p1ML=plot(MLt1{i},MLX1{i},'color',colors(i,:));

p1M=plot(Mt1{i},MX1{i},'color',colors(i,:));
% legend('Memoryless','Memory=0.2')
% set(get(get(p1R,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

p1R=plot(Rt1{i},RX1{i},'color',colors(i,:));
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
   p4=plot(Xall2{i},Q1{i},'color',colors(i,:));
hold on

set(p4,'LineWidth',3)
set(gca,'FontSize',14) 
end
xlabel('States')
ylabel('Potential energy')

subplot(2,2,4) % slope of potential landscape 
hold on
for i=1:2
   p4=plot(Xall2{i},Dx1{i},'color',colors(i,:));
hold on

set(p4,'LineWidth',3)
set(gca,'FontSize',14) 
end
xlabel('States')
ylabel('Growth rates')
yline(0,':')
axis([0 5 -1.5 .8 ])
