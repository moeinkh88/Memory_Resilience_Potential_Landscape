clc
clear 

global r K A B 

% coefficients
r=.8;
K=3;
A=.2; 
B=.6;

F=@funCP; % herbivory model equation function
JF=@JfunCP; % Jacobian of function

al=1:-.2:.8; % order of derivatives (memory = 0 and 0.2)
t0=0; % initial time
T=500; % final time

h=0.01; % computation step size

xmin=[0,1.9568]; % stable states of the model
xmax=.8432; % unstable state of the model

X0=xmin(2);

% perturbation
P=0.12; % P is strength of perturbation (B+p)
pt0=10; % the time when perturbation starts
pt1=20; % the time when perturbation ends

IndxP=pt1/h;

dt=xmin(2)-xmax; % distance to threshold from stable point
Eps=dt*1e-4; % scale for recovery rate

%% evaluating Resilience
for i=1:length(al)
        
[t,x] = FDE_PI2_IM(al(i),F,JF,t0,T,X0,h,[P,pt0,pt1]);% solving the model to find states
Xmin=x(pt0/h-1);
RX(i,:)=x;
ind=find(abs(x(IndxP:end)-Xmin)<Eps);

% EngRL(i)=t(indx(1))+pt1;
EngRL(i)=(2*(abs(Xmin-min(x)))/(abs(Xmin-min(x))+Eps)-1)/(t(ind(1)));
end
%% evaluating the potential landscape when there is no perturbation
% converting model to polynomial format
a1=A*r-B;
a2=r*(1-A/K);
a3=-r/K;
Dlt=a2^2-4*a1*a3;
r1=(-a2-sqrt(Dlt))/(2*a3);
r2=(-a2+sqrt(Dlt))/(2*a3);

par=[a3,a2,a1];
H=.0001;
xt=xmax-.1:H:xmax+xmin(2);
U=-polyint(par);
Q=polyval(U,xt); % potential function values
TF=islocalmin(Q);
TF1=islocalmax(Q);
Hp=Q(TF1); % height of potential landscape on threshold (peak)
Lp=Q(TF);  % height of potential landscape in attractions
dp=abs(Hp-Lp) % dept of basin of attractions (left and right)

DQ=diff(Q)/H;
indx=xt(TF);
XRattr=find(xt==indx);
Curv=-(DQ(XRattr+1)-DQ(XRattr))/H % finding the curvation of the stable state

%% plotting the dynamics
colors=[0 0.4470 0.7410
    0.8500 0.3250 0.0980;
0.4660 0.6740 0.1880
0.9290 0.6940 0.1250];

figure

 % Highlight background as perturbation
    vb = [pt0 min(min(RX(:,:))); pt1 min(min(RX(:,:))); pt1 max(max(RX(:,:))); pt0 max(max(RX(:,:)))];
    f = [1 2 3 4];
    h11=patch('Faces',f,'Vertices',vb,'FaceColor',colors(4,:),'EdgeColor','non', 'FaceAlpha',.16);
%     set(gca,'YScale','log')
    hold on


p=plot(t,RX(1,:),'color',colors(1,:));
hold on

set(p,'LineWidth',4)
set(gca,'FontSize',14)

p=plot(t,RX(2,:),'color',colors(2,:));

set(p,'LineWidth',3)
set(gca,'FontSize',14)

legend(['P=',num2str(P)],'Memory=0','Memory=0.2')
xlabel('Time')
ylabel('States')

axis([0 50 min(min(RX(:,:))) max(max(RX(:,:)))])


RelEff=diff(EngRL)/sum(EngRL)
DiffEngRLrate=abs(diff(EngRL))%difference of resilience

%% movement in potential landscape
par=[P,pt0,pt1];

[Q1LM,Xall1LM]=solveFDEM1(par,al(2)); % Potential with perturbation left side
[Q1RM,Xall1RM]=solveFDEM(par,al(2)); % Potential with perturbation right side
%% Evaluating potential during and after perturbation

B=B+P;
a1=A*r-B;
a2=r*(1-A/K);
a3=-r/K;
par=[a3,a2,a1];
Dlt=a2^2-4*a1*a3;
r1=(-a2-sqrt(Dlt))/(2*a3);
r2=(-a2+sqrt(Dlt))/(2*a3);

if r1>0 && r2>0 && isreal(r1) && isreal(r2)
    [Xall1,Q1]=PotenDepth(par,al(1)); % Potential without perturbation
else
    [Q1,Xall1]=solveFDE(par,al(1)); % Potential with perturbation
end

B=B-P;
a1=A*r-B;
a2=r*(1-A/K);
a3=-r/K;
par=[a3,a2,a1];
Dlt=a2^2-4*a1*a3;
r1=(-a2-sqrt(Dlt))/(2*a3);
r2=(-a2+sqrt(Dlt))/(2*a3);

[Xall,Q]=PotenDepth(par,al(1)); % Potential without perturbation
[XallM,QM]=PotenDepth(par,al(2)); % Potential without perturbation

%% plotting the potential before the perturbation


figure 
p1=plot(Xall,Q,'color',colors(1,:));
hold on

%shifting the plot to the level of fixed initial point
TF=islocalmin(Q); adj=QM(TF)-Q(TF);

p2=plot(XallM,QM-adj,'color',colors(2,:));
set(p1,'LineWidth',4)
set(p2,'LineWidth',4)
set(gca,'FontSize',14)

ylim([-.05 .12])
ylabel('Potential energy')
xlabel('States')
%% plotting the potential during the perturbation
figure

ind0=find(abs(Xall1-Xall(TF))<.006);ind0=ind0(1);
Q1=Q1-(Q1(ind0)-Q(TF));% adjusting fixed initial point

% p1fade=plot(Xall,Q,'color',colors(1,:),'LineWidth',4); p1fade.Color(4)=.2;
hold on
p1=plot(Xall1,Q1,'color',colors(1,:));
ind=find(abs(Xall1-Xall1LM(end))<.01);ind=ind(1);
QLM=Q1LM-(Q1LM(end)-Q1(ind));

% p1fade=plot(XallM,QM-adj,'color',colors(2,:),'LineWidth',4); p1fade.Color(4)=.2;
p2=plot(Xall1LM,QLM,'color',colors(2,:));
set(p1,'LineWidth',4)
set(p2,'LineWidth',4)
set(gca,'FontSize',14)

xlim([0 2.5])
ylim([-.3 .25])

ylabel('Potential energy')
xlabel('States')

%% plotting the potential after the perturbation
figure

ind3=find(abs(Xall1-1.43455)<.006);ind3=ind3(1);
p1fade=plot(Xall1(ind3:ind0),Q1(ind3:ind0),'color',colors(1,:),'LineWidth',4); p1fade.Color(4)=.2;

hold on

ind2=find(abs(Xall-1.43455)<.006);ind2=ind2(1);
QR=Q-(Q(ind2)-Q1(ind3));%adjustment
p1=plot(Xall,QR,'color',colors(1,:));

ind=find(abs(Xall1LM-1.57736)<.005);ind=ind(1);
p2fade=plot(Xall1LM(ind:end),QLM(ind:end),'color',colors(2,:),'LineWidth',4); p2fade.Color(4)=.2;

QRM=Q1RM-(Q1RM(end)-QLM(ind));%adjustment
p2=plot(Xall1RM,QRM,'color',colors(2,:));


set(p1,'LineWidth',4)
set(p2,'LineWidth',4)
set(gca,'FontSize',14)

ylabel('Potential energy')
xlabel('States')
% xlim([1 2.2])
ylim([-.16 .01])



