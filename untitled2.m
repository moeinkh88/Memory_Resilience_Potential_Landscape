clc
clear

load("CurvL.mat");
load("CurvR.mat");
load("EcoRL.mat")
load("EngRL.mat")
load("Hpp.mat")
load("Ldp.mat")
load("LLp.mat")
load("p.mat")
load("P08.mat")
load("PAR.mat")
load("Potential.mat")
load("Rdp.mat")
load("RLp.mat")
load("States.mat")
load("xmax.mat")
load("xmin.mat")

xt=0:.00001:5; % steps for polynomial function
Order=1:-.02:.8; % order of derivatives (1-memory) 

Il=1;Jl=1;Ir=1;Jr=1;

Nsample=1000;
for i=1:Nsample

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


[R,pv1]=corrcoef(Ldp(i,:),1-Order);
corrLdp(i)=R(2);
[R,pv2]=corrcoef(Rdp(i,:),1-Order);
corrRdp(i)=R(2);
[R,pv3]=corrcoef(LLp(i,:),1-Order);
corrLLp(i)=R(2);
[R,pv4]=corrcoef(RLp(i,:),1-Order);
corrRLp(i)=R(2);
[R,pv5]=corrcoef(Hpp(i,:),1-Order);
corrHpp(i)=R(2);
[R,pv6]=corrcoef(CurvL(i,:),1-Order);
corrCurvL(i)=R(2);
[R,pv7]=corrcoef(CurvR(i,:),1-Order);
corrCurvR(i)=R(2);


[RHO,pv8] = corr(CurvL(i,:)',1-Order','Type','Spearman');
corrCurvLSpear(i)=RHO;
[RHO,pv9] = corr(CurvR(i,:)',1-Order','Type','Spearman');
corrCurvRSpear(i)=RHO;

if corrRdp(i)>0
    pRelRdp(Ir)=abs(CurvR(i,1)-CurvR(i,end))/abs(CurvR(i,1));
    Ir=Ir+1;
else
    nRelRdp(Jr)=abs(CurvR(i,1)-CurvR(i,end))/abs(CurvR(i,1));
    Jr=Jr+1;
end

if corrLdp(i)>0
    pRelLdp(Il)=abs(CurvL(i,1)-CurvL(i,end))/abs(CurvL(i,1));
    Il=Il+1;
else
    nRelLdp(Jl)=abs(CurvL(i,1)-CurvL(i,end))/abs(CurvL(i,1));
    Jl=Jl+1;
end

pv(i)=mean([pv1(2),pv2(2),pv3(2),pv4(2),pv5(2),pv8,pv9]); %avarage of pvalues

toc
end
%%


Cpr=[0.4940 0.1840 0.5560];
Cgr=[0.4660 0.6740 0.1880];

figure

indPc=find(corrRdp>0);
indNc=find(corrRdp<0);

p3=scatter(indNc,corrRdp(indNc),[],'k');
hold on
p33=scatter(indPc,corrRdp(indPc),[],'k');
set(p3,'LineWidth',.5)
set(p33,'LineWidth',.5)
set(gca,'FontSize',14)
ylabel('Corr pd and memory')
xlabel('Sample models')

%violin plot
figure
indP1 = repmat({'Corr pd and memory>0'},length(indPc),1);
indP2 = repmat({'Corr pd and memory<0'},length(indNc),1);

violinplot([corrRdp(indPc)';corrRdp(indNc)'],[indP1;indP2])% 
ylabel('Corr pd and memory')
set(gca,'FontSize',14)

%% Effects of memory on curvature

figure

p2=scatter(1:Nsample,corrCurvR,[],'k');
% p2=plot(alpha0,corrLdp);
set(p2,'LineWidth',.5)
set(gca,'FontSize',14)
ylabel('Corr curv and memory (Pearson)')
xlabel('Sample models')

figure

p2=scatter(1:Nsample,corrCurvRSpear,[],'k');
% p2=plot(alpha0,corrLdp);
set(p2,'LineWidth',.5)
set(gca,'FontSize',14)
ylabel('Corr curv and memory (Spearman)')
xlabel('Sample models')

%violin plot
figure
indP11 = repmat({'Pearson'},length(corrCurvR),1);
indP22 = repmat({'Spearman'},length(corrCurvRSpear),1);

violinplot([corrCurvR';corrCurvRSpear'],[indP11;indP22])% 
ylabel('Corr curv and memory')
set(gca,'FontSize',14)

figure
g3 = repmat({'Memory=0 R valley'},length(CurvR(:,1)),1);
g4 = repmat({'Memory=0.2 R valley'},length(CurvR(:,end)),1);

violinplot([CurvR(:,1);CurvR(:,end)],[g3;g4]);
ax = gca;
ax.YAxis.Scale ="log";

ylabel('Curvature')


%% Memory effects on potential depth in terms of potential depth and curvature and initial slope of potential


[CRcurv,iaRcurv] = unique(CurvR(:,1));
figure

indPc=find(corrRdp(iaRcurv)>0);
indNc=find(corrRdp(iaRcurv)<0);

plot(CRcurv(indNc),corrRdp(iaRcurv(indNc)),'o','color',Cpr)
hold on
plot(CRcurv(indPc),corrRdp(iaRcurv(indPc)),'o','color',Cgr)
xlabel('Curvature')
ylabel('Corr pd and memory')

set(gca,'FontSize',14)

[CRdp,iaRdp] = unique(Rdp(:,1));
figure

indPc1=find(corrRdp(iaRdp)>0);
indNc1=find(corrRdp(iaRdp)<0);

semilogx(CRdp(indNc1),corrRdp(iaRdp(indNc1)),'o','color',Cpr)
hold on
semilogx(CRdp(indPc1),corrRdp(iaRdp(indPc1)),'o','color',Cgr)
xlabel('Potential depth')
ylabel('Corr pd and memory')

set(gca,'FontSize',14)
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

%% scatter curv

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
%% scatter DeltaDp
figure

scatter(Ldp(:,1)-Ldp(:,11),Rdp(:,1)-Rdp(:,11),'k')
hold on
scatter(Ldp(:,1)-Ldp(:,2),Rdp(:,1)-Rdp(:,2),'m')
set(gca,'xscale','log')
set(gca,'yscale','log')
% symlog(gca,'xy',-1.7)
% yl = get(gca,'ytick');
% set(gca,'yticklabel',sign(yl).*10.^abs(yl))
% xl = get(gca,'xtick');
% set(gca,'xticklabel',sign(xl).*10.^abs(xl))

%% delta curv vs dp

figure

hold on

indP=find((Rdp(:,1)-Rdp(:,11))>0);
indN=find((Rdp(:,1)-Rdp(:,11))<0);
scatter(Rdp(indP,1),abs((CurvR(indP,1)-CurvR(indP,11))./CurvR(indP,1)),[],Cpr)
scatter(Rdp(indN,1),abs((CurvR(indN,1)-CurvR(indN,11))./CurvR(indN,1)),[],Cgr)

set(gca,'xscale','log')
% set(gca,'yscale','log')
xlabel('Potential depth')
ylabel('Relative effect of memory on curvature')

set(gca,'FontSize',14)
%% delta curv vs curv

figure

hold on

indP=find((Rdp(:,1)-Rdp(:,11))>0);
indN=find((Rdp(:,1)-Rdp(:,11))<0);
scatter(CurvR(indP,1),abs((CurvR(indP,1)-CurvR(indP,11))./CurvR(indP,1)),[],Cpr)
scatter(CurvR(indN,1),abs((CurvR(indN,1)-CurvR(indN,11))./CurvR(indN,1)),[],Cgr)

% set(gca,'xscale','log')
% set(gca,'yscale','log')
xlabel('Curvature')
ylabel('Relative effect of memory on curvature')

set(gca,'FontSize',14)
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
%% memory effects on curvature in terms of + and - effects on potential depth
figure
g1 = repmat({'Corr pd and memory>0'},length(pRelRdp),1);
g2 = repmat({'Corr pd and memory<0'},length(nRelRdp),1);

boxplot([pRelRdp';nRelRdp'],[g1;g2])% 
ylabel('Relative effect of memory on curvature')

set(gca,'FontSize',14)

%% test

a3=PAR(:,1);
a2=PAR(:,2);
a1=PAR(:,3);

a2a1a3=a2.^2-4*a1.*a3;

[CRdisL,iadisL] = unique(a2a1a3);
plot(CRdisL,corrLdp(iadisL),'o')
xlabel('a3 (L valley)')
ylabel('Corr pd and memory')

%% test

a3=PAR(:,1);
a2=PAR(:,2);
a1=PAR(:,3);

indx1=find(a2>0);
indx2=find(a1(indx1)<0);

find(a2(indx2)>sqrt(4*a1(indx2).*a3(indx2)))
a2a1a3=a2.^2-4*a1.*a3;

[CRdisL,iadisL] = unique(a2a1a3);
plot(CRdisL,corrLdp(iadisL),'o')
xlabel('a3 (L valley)')
ylabel('Corr pd and memory')
