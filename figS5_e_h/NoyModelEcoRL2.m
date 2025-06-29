clc
clear 

global r K A B 

% coefficients
r=.8;
K=3;
A=.2;
B=.6;

p1=.187; % initial perturbation

F=@funCP;
JF=@JfunCP;

al=1:-.2:.8;
t0=0;
T=100;

h=0.01;

xmin=[0,1.9568];
xmax=.8432;

RX0=xmin(2);
pt0=10;pt1=20;

%%

for i=1:2
        tic
        P=p1; %perturb
while 1
[Rt,Rx] = FDE_PI2_IM(al(i),F,JF,t0,T,RX0,h,[P,pt0,pt1]);

Df=diff(Rx(:,end-1:end)');

if Rx(end)>xmax 
    P=P+0.0005;
    RX(i,:)=Rx;
else

    EcoRL(i)=(P-0.0005);  
%     EcoRL(i)=P;  
% RX(i,:)=Rx;

    break
end
toc
end
end

%%
colors=[0 0.4470 0.7410
    0.8500 0.3250 0.0980;
0.4660 0.6740 0.1880
0.9290 0.6940 0.1250];

figure
set(gcf,'renderer','Painters')

subplot(2,1,1)
 % Highlight background as perturbation
%     vb = [10 min(min(RX(:,:))); 20 min(min(RX(:,:))); 20 max(max(RX(:,:))); 10 max(max(RX(:,:)))];
        vb = [10 0; 20 0; 20 max(max(RX(:,:))); 10 max(max(RX(:,:)))];
    f = [1 2 3 4];
    h11=patch('Faces',f,'Vertices',vb,'FaceColor',colors(4,:),'EdgeColor','non', 'FaceAlpha',.16);
%     set(gca,'YScale','log')
    hold on


p=plot(Rt,RX(1,:),'color',colors(1,:));
hold on
yline(.8432,':')
text(5,.8432,'X_{u}','FontSize',14)
legend(['P=',num2str(EcoRL(1))])
title('Memory=0')

set(p,'LineWidth',3)
set(gca,'FontSize',14)


subplot(2,1,2)
 % Highlight background as perturbation
%     vb = [10 min(min(RX(:,:))); 20 min(min(RX(:,:))); 20 max(max(RX(:,:))); 10 max(max(RX(:,:)))];
        vb = [10 0; 20 0; 20 max(max(RX(:,:))); 10 max(max(RX(:,:)))];
    f = [1 2 3 4];
    h11=patch('Faces',f,'Vertices',vb,'FaceColor',colors(4,:),'EdgeColor','non', 'FaceAlpha',.46);
%     set(gca,'YScale','log')
    hold on


p=plot(Rt,RX(2,:),'color',colors(2,:));
hold on
yline(.8432,':')
text(5,.8432,'X_{u}','FontSize',14)
legend(['P=',num2str(EcoRL(2))])
title('Memory=0.2')

set(p,'LineWidth',3)
set(gca,'FontSize',14)
xlabel('Time')
ylabel('States')

RelEff=diff(EcoRL)/abs(EcoRL(1))
% DiffEngRLrate=diff(EcoRL)
% RelEff=diff(EcoRL)/diff(DP)
EcoRL

%%
par=[P-.0005,pt0,pt1];
% par=[P,pt0,pt1];

[Q1LM,Xall1LM]=solveFDEM1(par,al(2)); % Potential with perturbation left side
[Q1RM,Xall1RM]=solveFDEM(par,al(2)); % Potential with perturbation right side
%%

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
%%

figure 
p1=plot(Xall,Q,'color',colors(1,:));
hold on

%shifting the plot to the level of fixed initial point
TF=islocalmin(Q); adj=QM(TF)-Q(TF);

p2=plot(XallM,QM-adj,'color',colors(2,:));
set(p1,'LineWidth',4)
set(p2,'LineWidth',4)
set(gca,'FontSize',14)
%%
figure

ind0=find(abs(Xall1-Xall(TF))<.006);ind0=ind0(1);
Q1=Q1-(Q1(ind0)-Q(TF));% adjusting fixed initial point

% p1fade=plot(Xall,Q,'color',colors(1,:),'LineWidth',4); p1fade.Color(4)=.2;
hold on
p1=plot(Xall1,Q1,'color',colors(1,:));
ind=find(abs(Xall1-Xall1LM(end))<.005);ind=ind(1);
QLM=Q1LM-(Q1LM(end)-Q1(ind));

% p1fade=plot(XallM,QM-adj,'color',colors(2,:),'LineWidth',4); p1fade.Color(4)=.2;
p2=plot(Xall1LM,QLM,'color',colors(2,:));
set(p1,'LineWidth',4)
set(p2,'LineWidth',4)
set(gca,'FontSize',14)
%%
figure

ind3=find(abs(Xall1-.84291)<.006);ind3=ind3(1);
p1fade=plot(Xall1(ind3:ind0),Q1(ind3:ind0),'color',colors(1,:),'LineWidth',4); p1fade.Color(4)=.2;

hold on

ind2=find(abs(Xall-.84291)<.006);ind2=ind2(1);
QR=Q-(Q(ind2)-Q1(ind3));%adjustment
p1=plot(Xall,QR,'color',colors(1,:));

ind=find(abs(Xall1LM-.60125)<.005);ind=ind(1);
p2fade=plot(Xall1LM(ind:end),QLM(ind:end),'color',colors(2,:),'LineWidth',4); p2fade.Color(4)=.2;

QRM=Q1RM-(Q1RM(end)-QLM(ind));%adjustment
p2=plot(Xall1RM,QRM,'color',colors(2,:));


set(p1,'LineWidth',4)
set(p2,'LineWidth',4)
set(gca,'FontSize',14)
xlabel('States')