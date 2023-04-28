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
X0=.01;
MLX0=3.4953;
MX0=3.49531;
RX0=6;
F=@fun5;
JF=@Jfun;
T=2000;

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

DdX=cat(1,dx,flip(dxML),dxM,flip(dxR));
Xxall=cat(1,X(:,1:end-1)',flip(MLX(:,1:end-1)'),MX(:,1:end-1)',flip(RX(:,1:end-1)'));

NX=50000;

subplot(1,2,2)
for i=1:length(al)

[C,ia,ic] = unique(Xxall(:,i));
Xall(:,i)=linspace(X0,RX0, NX);
DX(:,i) = interp1(C,DdX(ia,i),Xall(:,i),'pchip');

Q(:,i) = -cumtrapz(Xall(:,i),DX(:,i));
p4=plot(Xall(:,i),Q(:,i),'color',colors(i,:));
% , 'LineStyle', 'none','Marker','.');
hold on

set(p4,'LineWidth',3)
set(gca,'FontSize',14)

P(:,i)=exp(-Q(:,i)/4);
end
yline(0)
xlabel('states (X)')

ylabel('Potential landscape')

%% potential slop

H=((C(end)-C(1))/NX);
DQ=diff(Q)/H;

figure

for i=1:length(al)
    
    [TF,p] = islocalmin(Q(:,i));
    indx=Xall(TF,i);
    XLattr=find(Xall(:,i)==indx(1));
    XRattr=find(Xall(:,i)==indx(2));
    
p4=plot(Xall(1:end-1,i),-DQ(:,i),'color',colors(i,:));
% , 'LineStyle', 'none','Marker','.');
hold on

set(p4,'LineWidth',3)
set(gca,'FontSize',14)

LamAttrL(i)=-(DQ(XLattr+1,i)-DQ(XLattr,i))/H;
LamAttrR(i)=-(DQ(XRattr+1,i)-DQ(XRattr,i))/H;
end
yline(0)
xlabel('States (X)')
axis tight
xline(indx(1))
xline(indx(2))
ylabel('Potential slop \newline rate of change')


text(1.5,.25,'recovery rate', 'FontSize',12)
text(4,.25,'recovery rate', 'FontSize',12)
text(1.5,.2,num2str(LamAttrL(1)), 'FontSize',12,'color',colors(1,:))
text(1.5,.15,num2str(LamAttrL(2)), 'FontSize',12,'color',colors(2,:))
text(4,.2,num2str(LamAttrR(1)), 'FontSize',12,'color',colors(1,:))
text(4,.15,num2str(LamAttrR(2)), 'FontSize',12,'color',colors(2,:))

%% curvature, analytical
% f(x,r,k,l,a)=(r.*x.*(1-x/k)-l*x.^2./(x.^2+a^2));
% curv=-f';
figure
xxx=X0:.01:RX0;
Curv=rr.*(xxx./KK-1)+rr/KK.*xxx +2.*lambda.*xxx./(xxx.^2+aa^2)-2.*lambda*xxx.^3./((xxx.^2+aa^2).^2);
 
plot(xxx,-Curv)
hold on
%% curvature

DDQ=-diff(DQ)/H;

% figure

for i=1:length(al)

p224=plot(Xall(1:end-2,i),DDQ(:,i),'color',colors(i,:));
% , 'LineStyle', 'none','Marker','.');
hold on

set(p224,'LineWidth',3)
set(gca,'FontSize',14)

end
xline(indx(1))
xline(indx(2))
yline(0)
xlabel('X(t)')
axis([0,6,-.2,.4])

ylabel('curvature')

%% mean exit time
figure

for i=1:length(al)
 TF1 = islocalmax(Q(:,i));
 [TF,p] = islocalmin(Q(:,i));
    indx=Xall(TF,i);
    indx1=Xall(TF1,i);
    Xmax=find(Xall(:,i)==indx1);
   
   MeanExit(:,i)=2*pi./sqrt(abs(DDQ(:,i).*DDQ(TF1,i))).*exp(abs(Q(Xmax,i)-Q(1:end-2,i))/0.1);
   p324=plot(Xall(1:end-2,i),MeanExit(:,i),'color',colors(i,:));

hold on

set(p324,'LineWidth',3)
set(gca,'FontSize',14)

xline(indx(1))
xline(indx(2))
    axis([0,6,0,400])
end
    
    
    
    
    
    
