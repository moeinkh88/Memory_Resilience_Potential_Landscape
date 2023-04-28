clc
clear
global rr KK lambda aa h T
rr=1;
KK=10;
lambda=2.75;
aa=1.6;
t0=0;
h=0.01;
al=1:-.2:.8;
X0=.5;
MLX0=3.49;
MX0=3.5;
RX0=6;
F=@fun5;
JF=@Jfun;
T=1000;

for i=1:length(al)
         [t,x] = FDE_PI2_IM(al(i),F,JF,t0,T,X0,h); 
%          [t,x] = FDE_PI12_PC(al(i),F,t0,T,X0,h); 
    X(i,:)=x;
end

for i=1:length(al)
         [MLt,MLx] = FDE_PI2_IM(al(i),F,JF,t0,T,MLX0,h);
%          [MLt,MLx] = FDE_PI12_PC(al(i),F,t0,T,MLX0,h); 
    MLX(i,:)=MLx;
end

for i=1:length(al)
         [Mt,Mx] = FDE_PI2_IM(al(i),F,JF,t0,T,MX0,h); 
%          [Mt,Mx] = FDE_PI12_PC(al(i),F,t0,T,MX0,h); 
    MX(i,:)=Mx;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Right
for i=1:length(al)
         [Rt,Rx] = FDE_PI2_IM(al(i),F,JF,t0,T,RX0,h); 
%          [Rt,Rx] = FDE_PI12_PC(al(i),F,t0,T,RX0,h); 
    RX(i,:)=Rx;
end

%%
colors=[0.4940 0.1840 0.5560
0.4660 0.6740 0.1880
0.6350 0.0780 0.1840];
figure
subplot(1,2,1)


for i=1:length(al)
p1=plot(t,X(i,:),'color',colors(i,:));
hold on


p1ML=plot(MLt,MLX(i,:),'color',colors(i,:));

p1M=plot(Mt,MX(i,:),'color',colors(i,:));
% legend('Memoryless','Memory=0.2')
% set(get(get(p1R,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

p1R=plot(Rt,RX(i,:),'color',colors(i,:));
legend('Memoryless','Memory=0.2','Memory=0.4')
set(get(get(p1R,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p1M,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p1ML,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


set(p1,'LineWidth',3)
set(gca,'FontSize',14)
set(p1ML,'LineWidth',3)
set(gca,'FontSize',14)
set(p1M,'LineWidth',3)
set(gca,'FontSize',14)
set(p1R,'LineWidth',3)
set(gca,'FontSize',14)
end

xlabel('time')
ylabel('X(t)')
% axis([0,100,X0,RX0])
%%potential landscape



dx=diff(X')./h;
dxML=diff(MLX')./h;
dxM=diff(MX')./h;
dxR=diff(RX')./h;

J=3:length(X)-2;
dx(J-1,:)=   1/(12*h).*(X(:,J-2)' - 8.*X(:,J-1)' + 8.*X(:,J+1)' - X(:,J+2)');
dxML(J-1,:)=   1/(12*h).*(MLX(:,J-2)' - 8.*MLX(:,J-1)' + 8.*MLX(:,J+1)' - MLX(:,J+2)');
dxM(J-1,:)=   1/(12*h).*(MX(:,J-2)' - 8.*MX(:,J-1)' + 8.*MX(:,J+1)' - MX(:,J+2)');
dxR(J-1,:)=   1/(12*h).*(RX(:,J-2)' - 8.*RX(:,J-1)' + 8.*RX(:,J+1)' - RX(:,J+2)');

DX=cat(1,dx,flip(dxML),dxM,flip(dxR));
Xall=cat(1,X(:,1:end-1)',flip(MLX(:,1:end-1)'),MX(:,1:end-1)',flip(RX(:,1:end-1)'));

subplot(1,2,2)
for i=1:length(al)


Q(:,i) = cumtrapz(Xall(:,i),DX(:,i));
p4=plot(Xall(:,i),-Q(:,i),'color',colors(i,:));
% , 'LineStyle', 'none','Marker','.');
hold on

set(p4,'LineWidth',3)
set(gca,'FontSize',14)

P(:,i)=exp(Q(:,i)/4);
end
yline(0)
xlabel('states (X)')

ylabel('Potential landscape')

%% probability

figure

PD1=P(:,1)/max(P(:,1));
PD2=P(:,2)/max(P(:,2));

pp4=plot(Xall(:,1), PD1,'color',colors(1,:));
hold on 
pp5=plot(Xall(:,2), PD2,'color',colors(2,:));

xline(MLX(2,1),'--')
xline(MX(2,1),'--')

set(pp4,'LineWidth',3)
set(gca,'FontSize',14)
set(pp5,'LineWidth',3)
set(gca,'FontSize',14)

xlabel('states (X)')

ylabel('Probability distribution P(X)')

NL=length(dx)+length(dxML);

SK1L=skewness(Xall(1:NL,1));
SK1R=skewness(Xall(NL+1:end,1));
SK2L=skewness(Xall(1:NL,2));
SK2R=skewness(Xall(NL+1:end,2));

var1L=var(Xall(1:NL,1));
var1R=var(Xall(NL+1:end,1));
var2L=var(Xall(1:NL,2));
var2R=var(Xall(NL+1:end,2));

md1L=median(Xall(1:NL,1),1);
md1R=median(Xall(NL+1:end,1),1);
md2L=median(Xall(1:NL,2),1);
md2R=median(Xall(NL+1:end,2),1);
% 

axis([0,6,.94,1])

text(3.7,.97,'Var', 'FontSize',12)
text(3.7,.96,'Sk', 'FontSize',12)
text(3.7,.95,'Med', 'FontSize',12)
text(4.2,.97,num2str(var1R), 'FontSize',12,'color',colors(1,:))
text(4.2,.96,num2str(SK1R), 'FontSize',12,'color',colors(1,:))
text(4.2,.95,num2str(md1R), 'FontSize',12,'color',colors(1,:))
text(4.2,.965,num2str(var2R), 'FontSize',12,'color',colors(2,:))
text(4.2,.955,num2str(SK2R), 'FontSize',12,'color',colors(2,:))
text(4.2,.945,num2str(md2R), 'FontSize',12,'color',colors(2,:))


text(.85,.97,'Var', 'FontSize',12)
text(.85,.96,'Sk', 'FontSize',12)
text(.85,.95,'Med', 'FontSize',12)
text(1.25,.97,num2str(var1L), 'FontSize',12,'color',colors(1,:))
text(1.25,.96,num2str(SK1L), 'FontSize',12,'color',colors(1,:))
text(1.3,.95,num2str(md1L), 'FontSize',12,'color',colors(1,:))
text(1.25,.965,num2str(var2L), 'FontSize',12,'color',colors(2,:))
text(1.25,.955,num2str(SK2L), 'FontSize',12,'color',colors(2,:))
text(1.3,.945,num2str(md2L), 'FontSize',12,'color',colors(2,:))

%% curvature of potential
% figure
% 
% XXX=cat(1,X(:,2000:end)',flip(MLX(:,1:end)'),MX(:,1:end)',flip(RX(:,1000:end)'));
% 
% for i=1:length(al)
% XattrL=MLX(i,end);
% XattrR=MX(i,end);
% 
% LamAttrL(1)=(dxML(end-1,i)-dxML(end,i))/h;
% LamAttrR(1)=(dxM(end-1,i)-dxM(end,i))/h;
% 
% % XL=cat(2,X(i,:),MLX(i,:));
% % XR=cat(2,MX(i,:),RX(i,:));
% tauReturnL1(i,:)=log(XattrL./X(i,2000:end))./LamAttrL;
% tauReturnL2(i,:)=log(XattrL./MLX(i,1:end))./LamAttrL;
% tauReturnR1(i,:)=log(XattrR./MX(i,1:end))./LamAttrR;
% tauReturnR2(i,:)=log(XattrR./RX(i,1000:end))./LamAttrR;
% 
% ReturnSpeed=cat(1,tauReturnL1(i,:)',flip(tauReturnL2(i,:)'),tauReturnR1(i,:)',flip(tauReturnR2(i,:)'));
% 
% end
% plot(XXX,ReturnSpeed)

%% potential slop

figure

for i=1:length(al)

p4=plot(Xall(:,i),DX(:,i),'color',colors(i,:));
% , 'LineStyle', 'none','Marker','.');
hold on

set(p4,'LineWidth',3)
set(gca,'FontSize',14)

Uprim1L=(dx(end-1,i)-dx(end,i))/h;
Uprim2L=(dxML(end-1,i)-dxML(end,i))/h;
LamAttrL(i)=mean([Uprim1L,Uprim2L]);
Uprim1R=(dxM(end-1,i)-dxM(end,i))/h;
Uprim2R=(dxR(end-1,i)-dxR(end,i))/h;
LamAttrR(i)=mean([Uprim1R,Uprim2R]);
end
yline(0)
xlabel('States (X)')
axis tight

ylabel(' Potential slop \newline rate of change')


text(1.5,.4,'recovery rate', 'FontSize',12)
text(4,.4,'recovery rate', 'FontSize',12)
text(1.5,.35,num2str(LamAttrL(1)), 'FontSize',12,'color',colors(1,:))
text(1.5,.3,num2str(LamAttrL(2)), 'FontSize',12,'color',colors(2,:))
text(4,.35,num2str(LamAttrR(1)), 'FontSize',12,'color',colors(1,:))
text(4,.3,num2str(LamAttrR(2)), 'FontSize',12,'color',colors(2,:))


%%

ddx=diff(dx)./h;
ddxML=diff(dxML)./h;
ddxM=diff(dxM)./h;
ddxR=diff(dxR)./h;

% J=3:length(dx)-2;
% ddx(J-1,:)=   1/(12*h).*(dx(J-2,:) - 8.*dx(J-1,:) + 8.*dx(J+1,:) - dx(J+2,:));
% ddxML(J-1,:)=   1/(12*h).*(dxML(J-2,:) - 8.*dxML(J-1,:) + 8.*dxML(J+1,:) - dxML(J+2,:));
% ddxM(J-1,:)=   1/(12*h).*(dxM(J-2,:) - 8.*dxM(J-1,:) + 8.*dxM(J+1,:) - dxM(J+2,:));
% ddxR(J-1,:)=   1/(12*h).*(dxR(J-2,:) - 8.*dxR(J-1,:) + 8.*dxR(J+1,:) - dxR(J+2,:));

DDX=cat(1,ddx,flip(ddxML),ddxM,flip(ddxR));
Xall2=cat(1,X(:,1:end-2)',flip(MLX(:,1:end-2)'),MX(:,1:end-2)',flip(RX(:,1:end-2)'));


figure

for i=1:length(al)

p224=plot(Xall2(:,i),DDX(:,i),'color',colors(i,:));
% , 'LineStyle', 'none','Marker','.');
hold on

set(p224,'LineWidth',3)
set(gca,'FontSize',14)

end

yline(0)
xlabel('X(t)')
axis tight

ylabel('curvature')


%% autocorrelation

