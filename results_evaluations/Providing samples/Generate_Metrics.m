% Creating samples among 1000 models with random parameters
clc
clear
%
global par
xt=0:.00001:5; % steps for polynomial function
Order=1:-.02:.8; % order of derivatives (1-memory) 

Il=1;Jl=1;Ir=1;Jr=1;

rng('default') % As of MATLAB R2014b and later, the default seed setting for rng('default') sets the seed to 0
s = rng;

Nsample=1000;
for i=1:Nsample
tic
    true=0;

while true==0
    a0=0;
    a1=-randi([0,4])*rand(1);
    a3=-randi([1,2])*rand(1);
    a2=sqrt(4*a1*a3)+randi([1,4])*rand(1);
par=[a3,a2,a1,a0]; % Randomly chosen coefficients of the polynomial     

Dlt=a2^2-4*a1*a3;

if a2>0 && a1<0 && a2>sqrt(Dlt)
   r1=(-a2-sqrt(Dlt))/(2*a3);
   r2=(-a2+sqrt(Dlt))/(2*a3);
    if r1<5 && r2<5
    
    U=-polyint(par); % potential of polynomial fun; -int(poly)
    U1=polyval(U,xt);
    TF=islocalmin(U1);
    TF1=islocalmax(U1);
    Depth(1)=abs(U1(TF1));
    Depth(2)=abs(U1(TF1)-U1(TF));
    if Depth(1)>.05 && Depth(2)>.05
    true=1;
    end
    end
end

end

PAR(i,:)=par;

[LRdp,Lp,Hp,LamAttr,Q,Xall]=solveFDE(par,Order,r1,r2);

Ldp(i,:)=LRdp(:,1);
Rdp(i,:)=LRdp(:,2);
LLp(i,:)=Lp(:,1);
RLp(i,:)=Lp(:,2);
Hpp(i,:)=Hp;
CurvL(i,:)=LamAttr(:,1);
CurvR(i,:)=LamAttr(:,2);
Potential{i,:,:}=Q;
States{i,:,:}=Xall;
%FLatness
Flat(i,:)=ceil(max(max(abs(CurvR))))+CurvR(i,:);
FlatL(i,:)=ceil(max(max(abs(CurvL))))+CurvL(i,:);
% find Xmin and Xmax
xmin(i,:)=xt(TF); % critical points: bottom of valley
xmax(i,:)=xt(TF1); % critical point: threshold 


% save('PAR.mat','PAR')
% save('Potential.mat','Potential')
% save('States.mat','States')
% save('Ldp.mat','Ldp')
% save('Rdp.mat','Rdp')
% save('LLp.mat','LLp')
% save('RLp.mat','RLp')
% save('Hpp.mat','Hpp')
% save('CurvL.mat','CurvL')
% save('CurvR.mat','CurvR')
% save('Flat.mat','Flat')
% save('FlatL.mat','FlatL')

toc
end

%% scatter dp

figure

hold on
p1=scatter(Ldp(:,1),Rdp(:,1));
p11=scatter(Ldp(:,end),Rdp(:,end),'+r');

set(gca,'xscale','log')
set(gca,'yscale','log')
legend('Memory=0','Memory=0.2')

set(p1,'LineWidth',1)
set(gca,'FontSize',14)
set(p11,'LineWidth',1)
set(gca,'FontSize',14)
ylabel('Potential depth of right valley')
xlabel('Potential depth of left valley')

%% scatter curvature

figure

hold on
p1=scatter(CurvL(:,1),CurvR(:,1));
p11=scatter(CurvL(:,end),CurvR(:,end),'+r');

set(gca,'xscale','log')
set(gca,'yscale','log')
legend('Memory=0','Memory=0.2')

% set(p1,'LineWidth',1)
set(gca,'FontSize',14)
% set(p11,'LineWidth',1)
set(gca,'FontSize',14)
ylabel('Curvatures on the right valley')
xlabel('Curvatures on the left valley')

%% dp vs curve

figure

hold on

p1=loglog(Rdp(:,1),CurvR(:,1),'ob');
% loglog(Ldp(:,1),CurvL(:,1),'ob')


p2=loglog(Rdp(:,11),CurvR(:,11),'+r');
% loglog(Ldp(:,11),CurvL(:,11),'+r')


% set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('Memory=0','Memory=0.2')

% yl = get(gca,'ytick');
% set(gca,'yticklabel',sign(yl).*10.^abs(yl))
set(gca,'xscale','log')
% set(gca,'yscale','log')
ylabel('Curvature')
xlabel('Potential depth')

set(p1,'LineWidth',.1)
set(p2,'LineWidth',.1)
set(gca,'FontSize',14)

axis tight

%% relative curve vs dp

figure

hold on

p1=loglog(Rdp(:,1),CurvR(:,1),'ob');
loglog(Ldp(:,1),CurvL(:,1),'ob')


p2=loglog(Rdp(:,11),CurvR(:,11),'+r');
loglog(Ldp(:,11),CurvL(:,11),'+r')


set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('Memory=0','Memory=0.2')

% yl = get(gca,'ytick');
% set(gca,'yticklabel',sign(yl).*10.^abs(yl))
set(gca,'xscale','log')
% set(gca,'yscale','log')
ylabel('Curvature')
xlabel('Potential depth')

axis tight
