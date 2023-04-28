% correlation between memory and potential landscape for polynomial
% function of degree 3
clc
clear
%
global par
xt=-10:.00001:10; % steps for polynomial function
Order=1:-.02:.8; % order of derivatives (1-memory) 

Il=1;Jl=1;Ir=1;Jr=1;
load('PAR')

Nsample=1000;
for i=1:length(PAR(:,1))
tic
par=PAR(i,:);
U=-polyint(par); % potential of polynomial fun; -int(poly)
U1=polyval(U,xt);
TF=islocalmin(U1);
TF1=islocalmax(U1);
xmin=xt(TF); % critical points: bottom of valley
xmax=xt(TF1); % critical point: threshold 
X0=xmin(1)-(xmax-xmin(1));
RX0=xmin(2)+(xmin(2)-xmax);
xtt=X0:.001:RX0;
U2=polyval(U,xtt);
DUL(i)=(U2(2)-U2(1))/.001;
DUR(i)=(U2(end)-U2(end-1))/.001;

[LRdp,Lp,Hp,LamAttr,Q,Xall]=solveFDE(U1,TF,TF1,Order,xt);

Ldp(i,:)=LRdp(:,1);
Rdp(i,:)=LRdp(:,2);
LLp(i,:)=Lp(:,1);
RLp(i,:)=Lp(:,2);
Hpp(i,:)=Hp;
CurvL(i,:)=LamAttr(:,1);
CurvR(i,:)=LamAttr(:,2);
Potential{i,:,:}=Q;
States{i,:,:}=Xall;

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

%% memory effects on curvature in terms of + and - effects on potential depth
figure
g1 = repmat({'+Corr pd, Left valley'},length(pRelLdp),1);
g2 = repmat({'-Corr pd, Left valley'},length(nRelLdp),1);
g3 = repmat({'+Corr pd, Right valley'},length(pRelRdp),1);
g4 = repmat({'-Corr pd, Right valley'},length(nRelRdp),1);

boxplot([pRelLdp';nRelLdp';pRelRdp';nRelRdp'],[g1;g2;g3;g4])
ylabel('Relative effect of memory on curvature')
%% memory effects in Curvature vs potential depth

figure

CurvREL=abs(CurvL(:,1)-CurvL(:,end))./abs(CurvL(:,1));
CurvRELR=abs(CurvR(:,1)-CurvR(:,end))./abs(CurvR(:,1));
RpdRELR=abs(Rdp(:,1)-Rdp(:,end))./abs(Rdp(:,1));
LpdRELR=abs(Ldp(:,1)-Ldp(:,end))./abs(Ldp(:,1));

scatter(LpdRELR,CurvREL)
hold on
scatter(RpdRELR,CurvRELR,'r')

ylabel('Relcurv')
xlabel('Reldp')

%% memory effects in Curvature vs potential depth in terms of + & - effects on potential depth
figure

for i=1:Nsample

if corrLdp(i)<0 && corrRdp(i)<0
    subplot(2,2,1)
    hold on
CurvREL=abs(CurvL(i,1)-CurvL(i,end))./abs(CurvL(i,1));
CurvRELR=abs(CurvR(i,1)-CurvR(i,end))./abs(CurvR(i,1));
RpdRELR=abs(Rdp(i,1)-Rdp(i,end))./abs(Rdp(i,1));
LpdRELR=abs(Ldp(i,1)-Ldp(i,end))./abs(Ldp(i,1));

plot(LpdRELR,CurvREL,'ob')
plot(RpdRELR,CurvRELR,'or')

ylabel('Relcurv')
xlabel('Reldp (-/-)')

    elseif corrLdp(i)<0 && corrRdp(i)>0
        subplot(2,2,2)
        hold on
        CurvREL=abs(CurvL(i,1)-CurvL(i,end))./abs(CurvL(i,1));
CurvRELR=abs(CurvR(i,1)-CurvR(i,end))./abs(CurvR(i,1));
RpdRELR=abs(Rdp(i,1)-Rdp(i,end))./abs(Rdp(i,1));
LpdRELR=abs(Ldp(i,1)-Ldp(i,end))./abs(Ldp(i,1));
        plot(LpdRELR,CurvREL,'ob')
plot(RpdRELR,CurvRELR,'or')

ylabel('Relcurv')
xlabel('Reldp (-/+)')

    elseif corrLdp(i)>0 && corrRdp(i)>0
    subplot(2,2,3)
    hold on
    CurvREL=abs(CurvL(i,1)-CurvL(i,end))./abs(CurvL(i,1));
CurvRELR=abs(CurvR(i,1)-CurvR(i,end))./abs(CurvR(i,1));
RpdRELR=abs(Rdp(i,1)-Rdp(i,end))./abs(Rdp(i,1));
LpdRELR=abs(Ldp(i,1)-Ldp(i,end))./abs(Ldp(i,1));
    plot(LpdRELR,CurvREL,'ob')
plot(RpdRELR,CurvRELR,'or')

ylabel('Relcurv')
xlabel('Reldp (+/+)')

    elseif corrLdp(i)>0 && corrRdp(i)<0
        subplot(2,2,4)
        hold on
        CurvREL=abs(CurvL(i,1)-CurvL(i,end))./abs(CurvL(i,1));
CurvRELR=abs(CurvR(i,1)-CurvR(i,end))./abs(CurvR(i,1));
RpdRELR=abs(Rdp(i,1)-Rdp(i,end))./abs(Rdp(i,1));
LpdRELR=abs(Ldp(i,1)-Ldp(i,end))./abs(Ldp(i,1));
        
plot(LpdRELR,CurvREL,'ob')
plot(RpdRELR,CurvRELR,'or')

ylabel('Relcurv')
xlabel('Reldp (+/-)')
end

end
title('blue=left, red=right')

%%  memory effects in Curvature vs potential depth in terms of + & - effects on potential depth
figure

for i=1:Nsample


if corrLLp(i)<0 && corrRLp(i)<0
    subplot(1,2,1)
    hold on
CurvREL=abs(CurvL(i,1)-CurvL(i,end))./abs(CurvL(i,1));
CurvRELR=abs(CurvR(i,1)-CurvR(i,end))./abs(CurvR(i,1));
RpdRELR=abs(Rdp(i,1)-Rdp(i,end))./abs(Rdp(i,1));
LpdRELR=abs(Ldp(i,1)-Ldp(i,end))./abs(Ldp(i,1));

plot(LpdRELR,CurvREL,'ob')
plot(RpdRELR,CurvRELR,'or')

ylabel('Relcurv')
xlabel('Reldp (Height Left/Right attr -/-)')

    elseif corrLLp(i)<0 && corrRLp(i)>0
        subplot(1,2,2)
        hold on
        CurvREL=abs(CurvL(i,1)-CurvL(i,end))./abs(CurvL(i,1));
CurvRELR=abs(CurvR(i,1)-CurvR(i,end))./abs(CurvR(i,1));
RpdRELR=abs(Rdp(i,1)-Rdp(i,end))./abs(Rdp(i,1));
LpdRELR=abs(Ldp(i,1)-Ldp(i,end))./abs(Ldp(i,1));
  
plot(LpdRELR,CurvREL,'ob')
plot(RpdRELR,CurvRELR,'or')

ylabel('Relcurv')
xlabel('Reldp (Height Left/Right attr -/+)')

end

end
title('blue=left, red=right')
%% correlation of memory and potential depth and Height of critical points

figure

subplot(2,3,2)
p8=scatter(1:Nsample,pv);
set(p8,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Average of all Pvalues')
xlabel('Sample models')

subplot(2,3,1)
p2=scatter(1:Nsample,corrLdp);
% p2=plot(alpha0,corrLdp);
set(p2,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr pd and memory (left valley)')
xlabel('Sample models')

subplot(2,3,3)
p3=scatter(1:Nsample,corrRdp);
% p3=plot(alpha0,corrLdp);
set(p3,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr pd and memory (right valley)')
xlabel('Sample models')

subplot(2,3,4)
p4=scatter(1:Nsample,corrLLp);
% p4=plot(alpha0,corrLdp);
set(p4,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr height of L attr and memory')
xlabel('Sample models')

subplot(2,3,6)
p5=scatter(1:Nsample,corrRLp);
% p5=plot(alpha0,corrLdp);
set(p5,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr height of R attr and memory')
xlabel('Sample models')

subplot(2,3,5)
p6=scatter(1:Nsample,corrHpp);
% p6=plot(alpha0,corrLdp);
set(p6,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr height of threshold and memory')
xlabel('Sample models')


%% Effects of memory on curvature

figure

subplot(3,2,1)
p2=scatter(1:Nsample,corrCurvL);
% p2=plot(alpha0,corrLdp);
set(p2,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr L curv and \newline memory (Pear)')
xlabel('Sample models')

subplot(3,2,2)
p2=scatter(1:Nsample,corrCurvR);
% p2=plot(alpha0,corrLdp);
set(p2,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr R curv and \newline memory (Pear)')
xlabel('Sample models')

subplot(3,2,3)
p2=scatter(1:Nsample,corrCurvLSpear);
% p2=plot(alpha0,corrLdp);
set(p2,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr L curv and \newline memory (Sp)')
xlabel('Sample models')

subplot(3,2,4)
p2=scatter(1:Nsample,corrCurvRSpear);
% p2=plot(alpha0,corrLdp);
set(p2,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr R curv and \newline memory (Sp)')
xlabel('Sample models')

subplot(3,2,[5,6])
g1 = repmat({'Memory=0 L valley'},length(CurvL(:,1)),1);
g2 = repmat({'Memory=0.2 L valley'},length(CurvL(:,end)),1);
g3 = repmat({'Memory=0 R valley'},length(CurvR(:,1)),1);
g4 = repmat({'Memory=0.2 R valley'},length(CurvR(:,end)),1);

boxplot([CurvL(:,1);CurvL(:,end);CurvR(:,1);CurvR(:,end)],[g1;g2;g3;g4]);
ax = gca;
ax.YAxis.Scale ="log";

ylabel('Curvature')
% scatter(1:Nsample,CurvL(:,1),'k');
% hold on
% scatter(1:Nsample,CurvL(:,end),'m');
% ylabel('left Curv (black NoMemory)')
% xlabel('samples')
% 
% subplot(2,2,4)
% scatter(1:Nsample,CurvR(:,1),'k');
% hold on
% scatter(1:Nsample,CurvR(:,end),'m');
% ylabel('right Curv (black NoMemory)')
% xlabel('samples')

%% Memory effects on potential depth in terms of potential depth and curvature and initial slope of potential

[CRU,iaRU] = unique(DUR);
[CLU,iaLU] = unique(DUL);
figure
subplot (3,2,2)
semilogx(CRU,corrRdp(iaRU),'o')
xlabel('Slope of potential (R valley)')
ylabel('Corr pd and memory')
subplot (3,2,1)
semilogx(CLU,corrLdp(iaLU),'o')
xlabel('Slope of potential (L valley)')
ylabel('Corr pd and memory')

[CRcurv,iaRcurv] = unique(CurvR(:,1));
[CLcurv,iaLcurv] = unique(CurvL(:,1));
subplot (3,2,4)
plot(CRcurv,corrRdp(iaRcurv),'o')
xlabel('Curv (R valley)')
ylabel('Corr pd and memory')
subplot (3,2,3)
plot(CLcurv,corrLdp(iaLcurv),'o')
xlabel('Curv (L valley)')
ylabel('Corr pd and memory')

[CRdp,iaRdp] = unique(Rdp(:,1));
[CLdp,iaLdp] = unique(Ldp(:,1));
subplot (3,2,6)
plot(CRdp,corrRdp(iaRdp),'o')
xlabel('Pd (R valley)')
ylabel('Corr pd and memory')
subplot (3,2,5)
plot(CLdp,corrLdp(iaLdp),'o')
xlabel('Pd (L valley)')
ylabel('Corr pd and memory')
%% scatter dp

figure

hold on
p1=scatter(Ldp(:,1),Rdp(:,1));
p11=scatter(Ldp(:,11),Rdp(:,end),'sr');

set(gca,'xscale','log')
set(gca,'yscale','log')
legend('Memory=0','Memory=0.2')

set(p1,'LineWidth',2)
set(gca,'FontSize',14)
set(p11,'LineWidth',2)
set(gca,'FontSize',14)
ylabel('Potential depth of right valley')
xlabel('Potential depth of left valley')
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

%% delta dp vs delta curv

figure

hold on

scatter(abs(Rdp(:,1)-Rdp(:,11)),abs((CurvR(:,1)-CurvR(:,11))./CurvR(:,1)),'b')
scatter((Rdp(:,1)-Rdp(:,11)),abs((CurvR(:,1)-CurvR(:,11))./CurvR(:,1)),'m')

scatter(abs(Ldp(:,1)-Ldp(:,11)),abs((CurvL(:,1)-CurvL(:,11))./CurvL(:,1)),'b')
scatter((Ldp(:,1)-Ldp(:,11)),abs((CurvL(:,1)-CurvL(:,11))./CurvL(:,1)),'m')
set(gca,'xscale','log')
set(gca,'yscale','log')
%% dp vs curve

figure

hold on

scatter(Rdp(:,1),CurvR(:,1),'b')
scatter(Ldp(:,1),CurvL(:,1),'b')

scatter(Rdp(:,11),CurvR(:,11),'r')
scatter(Ldp(:,11),CurvL(:,11),'r')

% yl = get(gca,'ytick');
% set(gca,'yticklabel',sign(yl).*10.^abs(yl))
set(gca,'xscale','log')
%%set(gca,'yscale','log')

%% Delta curve vs dp

figure

hold on

scatter(Rdp(:,1),abs(CurvR(:,11)-CurvR(:,1))./-CurvR(:,1),'r')
scatter(Ldp(:,1),abs(CurvL(:,11)-CurvL(:,1))./-CurvL(:,1),'b')

% yl = get(gca,'ytick');
% set(gca,'yticklabel',sign(yl).*10.^abs(yl))
% set(gca,'xscale','log')
% set(gca,'yscale','log')

%% memory effects on curvature in terms of + and - effects on potential depth
figure
g1 = repmat({'+Corr pd'},2*length(pRelLdp),1);
g2 = repmat({'-Corr pd'},2*length(nRelLdp),1);

boxplot([[pRelLdp';pRelRdp'];[nRelLdp';nRelRdp']],[g1;g2])% this is rel dp!!!!
ylabel('Relative effect of memory on curvature')


%% distinguish
for i=1:Nsample
tan=polyval(polyder(PAR(i,:)),roots(PAR(i,:)));

ipos=find(tan>0);

SlopThreShold(i)=tan(ipos);

end
figure
subplot (1,2,2)
[CRslop,iaslop] = unique(SlopThreShold);
semilogx(CRslop,corrRdp(iaslop),'o')
xlabel('slop tipping point (R valley)')
ylabel('Corr pd and memory')
subplot (1,2,1)
semilogx(CRslop,corrLdp(iaslop),'o')
xlabel('slop tipping point (L valley)')
ylabel('Corr pd and memory')

%% distinguish
for i=1:Nsample
rrr=roots(PAR(i,:));

DistAttr2ThreShold(i,:)=diff(unique(rrr));

end
figure
subplot (1,2,1)
[CRdisL,iadisL] = unique(DistAttr2ThreShold(:,1));
semilogx(CRdisL,corrLdp(iadisL),'o')
xlabel('DistAttr2ThreShold (L valley)')
ylabel('Corr pd and memory')
subplot (1,2,2)
[CRdisR,iadisR] = unique(DistAttr2ThreShold(:,2));
semilogx(CRdisR,corrRdp(iadisR),'o')
xlabel('DistAttr2ThreShold (R valley)')
ylabel('Corr pd and memory')

%% distinguish
for i=1:Nsample
HeightLR=polyval(PAR(i,:),unique(roots(polyder(PAR(i,:)))));
diffH(i)=diff(HeightLR);
end
figure
[CRdisL,iadisL] = unique(diffH);
semilogx(CRdisL,corrLdp(iadisL),'o')
xlabel('(L valley)')
ylabel('Corr pd and memory')

%% distinguish NOOO
for i=1:Nsample
a3(i)=PAR(i,1);
a2(i)=PAR(i,2);
a1(i)=PAR(i,3);
a0(i)=PAR(i,4);

end
figure
subplot (1,2,1)
[CRdisL,iadisL] = unique(a3);
plot(CRdisL,corrLdp(iadisL),'o')
xlabel('a3 (L valley)')
ylabel('Corr pd and memory')
subplot (1,2,2)
[CRdisR,iadisR] = unique(a3);
plot(CRdisR,corrRdp(iadisR),'o')
xlabel('a3 (R valley)')
ylabel('Corr pd and memory')

figure
subplot (1,2,1)
[CRdisL,iadisL] = unique(a2);
plot(CRdisL,corrLdp(iadisL),'o')
xlabel('a2 (L valley)')
ylabel('Corr pd and memory')
subplot (1,2,2)
[CRdisR,iadisR] = unique(a3);
plot(CRdisR,corrRdp(iadisR),'o')
xlabel('a2 (R valley)')
ylabel('Corr pd and memory')

figure
subplot (1,2,1)
[CRdisL,iadisL] = unique(a1);
plot(CRdisL,corrLdp(iadisL),'o')
xlabel('a1 (L valley)')
ylabel('Corr pd and memory')
subplot (1,2,2)
[CRdisR,iadisR] = unique(a1);
plot(CRdisR,corrRdp(iadisR),'o')
xlabel('a1 (R valley)')
ylabel('Corr pd and memory')

figure
subplot (1,2,1)
[CRdisL,iadisL] = unique(a0);
plot(CRdisL,corrLdp(iadisL),'o')
xlabel('a0 (L valley)')
ylabel('Corr pd and memory')
subplot (1,2,2)
[CRdisR,iadisR] = unique(a0);
plot(CRdisR,corrRdp(iadisR),'o')
xlabel('a0 (R valley)')
ylabel('Corr pd and memory')

%%
