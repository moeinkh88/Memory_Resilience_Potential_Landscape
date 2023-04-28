% correlation between memory and potential landscape for polynomial
% function of degree 3
clc
clear
%
global par
xt=-10:.001:10; % steps for polynomial function
Order=1:-.02:.8; % order of derivatives (1-memory) 

% alpha1=1.1:.1:1.6;
% Na1=length(alpha1);
% par=[-1 1 alpha1(1) -1]; % coefficients of the polynomial 

alpha1=-2:.1:-1.5;
Na1=length(alpha1);
par=[-1 1 3 alpha1(1)]; % coefficients of the polynomial 

figure

colors=[0.4940 0.1840 0.5560
0.4660 0.6740 0.1880
0.6350 0.0780 0.1840
0 0.4470 0.7410
0.8500 0.3250 0.0980
0.9290 0.6940 0.1250
0.3010 0.7450 0.9330
0.6350 0.0780 0.1840
0 1 0
1 0 0
0 0 0
];

for i=1:Na1
    
% par(3)=alpha1(i);
par(4)=alpha1(i);

U=-polyint(par); % potential of poly fun; -int(poly)

% p1=plot(xt,polyval(U,xt));

[LRdp,Lp,Hp,X1,U1,XN,UN,LamAttr]=solveFDE(U,Order,xt);
Ldp(i,:)=LRdp(:,1);
Rdp(i,:)=LRdp(:,2);
LLp(i,:)=Lp(:,1);
RLp(i,:)=Lp(:,2);
Hpp(i,:)=Hp;
CurvL(i,:)=LamAttr(:,1);
CurvR(i,:)=LamAttr(:,2);

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

pv(i)=mean([pv1(2),pv2(2),pv3(2),pv4(2),pv5(2),pv6(2),pv7(2)]); %avarage of pvalues

p1=plot(X1,U1,'color',colors(i,:));
set(p1,'LineWidth',3)
set(gca,'FontSize',14)

hold on
p22=plot(XN,UN,'--','color',colors(i,:));
set(p22,'LineWidth',3)
set(gca,'FontSize',14)
end

%%

ylabel('Potential landscape')
xlabel('States (X)')
axis('tight')

figure

subplot(2,3,1)
p2=plot(alpha1,corrLdp);
% p2=plot(alpha0,corrLdp);
set(p2,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr left dp')
xlabel('Coeficient of X^1')

subplot(2,3,3)
p3=plot(alpha1,corrRdp);
% p3=plot(alpha0,corrLdp);
set(p3,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr right dp')
xlabel('Coeficient of X^1')

subplot(2,3,4)
p4=plot(alpha1,corrLLp);
% p4=plot(alpha0,corrLdp);
set(p4,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr height of left attr')
xlabel('Coeficient of X^1')

subplot(2,3,6)
p5=plot(alpha1,corrRLp);
% p5=plot(alpha0,corrLdp);
set(p5,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr height of right attr')
xlabel('Coeficient of X^1')

subplot(2,3,5)
p6=plot(alpha1,corrHpp);
% p6=plot(alpha0,corrLdp);
set(p6,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr height of threshold')
xlabel('Coeficient of X^1')
%%
figure
plot(alpha1,pv)

figure
subplot(1,2,2)
p9=plot(corrCurvR);
set(p9,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('corr Curvatur right valley')
xlabel('Coeficient of X^1')

subplot(1,2,1)
p7=plot(corrCurvL);
set(p7,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('corr Curvatur left valley')
xlabel('Coeficient of X^1')