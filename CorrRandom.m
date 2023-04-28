% correlation between memory and potential landscape for polynomial
% function of degree 3
clc
clear
%
global par
xt=-30:.001:30; % steps for polynomial function
Order=1:-.02:.8; % order of derivatives (1-memory) 

figure

for i=1:30
true=0;

while true==0
par=[-2*rand(1),rand(1),2*rand(1),-2*rand(1)]; % Randomly chosen coefficients of the polynomial     
U=-polyint(par); % potential of polynomial fun; -int(poly)
U1=polyval(U,xt);
TF=islocalmin(U1);
TF1=islocalmax(U1);
Depth=abs(U1(TF1)-U1(TF));

if length(Depth)==2  
    if Depth(1)>.2 && Depth(2)>.2
%     p1=plot(xt,polyval(U,xt));
    true=1;
    end
end
end

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

[R,pv8]=corrcoef(CurvL(i,:),Ldp(i,:));
corrCrvLpd(i)=R(2);
[R,pv9]=corrcoef(CurvR(i,:),Rdp(i,:));
corrCrvRpd(i)=R(2);

if i==1
LDP=cat(1,Ldp(i,:));
RDP=cat(1,Rdp(i,:));
CRVL=cat(1,CurvL(i,:));
CRVR=cat(1,CurvR(i,:));
else
LDP=cat(2,LDP,Ldp(i,:));
RDP=cat(2,RDP,Rdp(i,:));
CRVL=cat(2,CRVL,CurvL(i,:));
CRVR=cat(2,CRVR,CurvR(i,:));
end
pv(i)=mean([pv1(2),pv2(2),pv3(2),pv4(2),pv5(2),pv6(2),pv7(2),pv8(2),pv9(2)]); %avarage of pvalues

subplot(2,1,1)
hold on
p1=plot(X1,U1);
set(p1,'LineWidth',3)
set(gca,'FontSize',14)

subplot(2,1,2)
hold on
p22=plot(XN,UN);
set(p22,'LineWidth',3)
set(gca,'FontSize',14)
end

%%



ylabel('Potential landscape')
xlabel('States (X)')
axis('tight')


figure

subplot(3,3,1)
p2=plot(corrLdp);
set(p2,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr left dp')
xlabel('Coeficient of X^1')

subplot(3,3,3)
p3=plot(corrRdp);
set(p3,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr right dp')
xlabel('Coeficient of X^1')

subplot(3,3,4)
p4=plot(corrLLp);
set(p4,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr hight of left attr')
xlabel('Coeficient of X^1')

subplot(3,3,6)
p5=plot(corrRLp);
set(p5,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr hight of right attr')
xlabel('Coeficient of X^1')

subplot(3,3,5)
p6=plot(corrHpp);
set(p6,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Corr hight of threshold')
xlabel('Coeficient of X^1')


subplot(3,3,7)
p7=plot(corrCurvL);
set(p7,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('corr Curvatur left valley')
xlabel('Coeficient of X^1')

subplot(3,3,9)
p9=plot(corrCurvR);
set(p9,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('corr Curvatur right valley')
xlabel('Coeficient of X^1')

subplot(3,3,8)
p8=plot(pv);
set(p8,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Average of Pvalues')
xlabel('Coeficient of X^1')
%%
figure

subplot(2,1,1)
% p21=plot(corrCrvLpd);
p21=scatter(CRVL,LDP);
set(p21,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Dp L')
xlabel('curv L')

subplot(2,1,2)
% p21=plot(corrCrvRpd);
p21=scatter(CRVR,RDP);
set(p21,'LineWidth',3)
set(gca,'FontSize',14)
ylabel('Dp R')
xlabel('curv R')