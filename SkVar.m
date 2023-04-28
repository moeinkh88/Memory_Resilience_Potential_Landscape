clc
clear
global rr KK lambda aa h T
rr=1;
KK=10;
lambda=2.75;
aa=1.6;
t0=0;
h=0.01;
al=1:-.02:.8;
X0=.001;
MLX0=3.4953;
MX0=3.49531;
RX0=6;
F=@fun5;
JF=@Jfun;
T=1000;

for i=1:length(al)
%          [t,x] = FDE_PI2_IM(al(i),F,JF,t0,T,X0,h); 
         [t,x] = FDE_PI12_PC(al(i),F,t0,T,X0,h); 
    X(i,:)=x;
end

for i=1:length(al)
%          [MLt,MLx] = FDE_PI2_IM(al(i),F,JF,t0,T,MLX0,h);
         [MLt,MLx] = FDE_PI12_PC(al(i),F,t0,T,MLX0,h); 
    MLX(i,:)=MLx;
end

for i=1:length(al)
%          [Mt,Mx] = FDE_PI2_IM(al(i),F,JF,t0,T,MX0,h); 
         [Mt,Mx] = FDE_PI12_PC(al(i),F,t0,T,MX0,h); 
    MX(i,:)=Mx;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Right
for i=1:length(al)
%          [Rt,Rx] = FDE_PI2_IM(al(i),F,JF,t0,T,RX0,h); 
         [Rt,Rx] = FDE_PI12_PC(al(i),F,t0,T,RX0,h); 
    RX(i,:)=Rx;
end

%%
colors=[0.4940 0.1840 0.5560
0.4660 0.6740 0.1880
0.6350 0.0780 0.1840];
% figure
% subplot(1,2,1)

% 
% for i=1:length(al)
% p1=plot(t,X(i,:),'color',colors(i,:));
% hold on
% 
% 
% p1ML=plot(MLt,MLX(i,:),'color',colors(i,:));
% 
% p1M=plot(Mt,MX(i,:),'color',colors(i,:));
% % legend('Memoryless','Memory=0.2')
% % set(get(get(p1R,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% 
% p1R=plot(Rt,RX(i,:),'color',colors(i,:));
% legend('Memoryless','Memory=0.2','Memory=0.4')
% set(get(get(p1R,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(p1M,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(p1ML,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% 
% 
% set(p1,'LineWidth',3)
% set(gca,'FontSize',14)
% set(p1ML,'LineWidth',3)
% set(gca,'FontSize',14)
% set(p1M,'LineWidth',3)
% set(gca,'FontSize',14)
% set(p1R,'LineWidth',3)
% set(gca,'FontSize',14)
% end
% 
% xlabel('time')
% ylabel('X(t)')



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

% subplot(1,2,2)
for i=1:length(al)


Q(:,i) = cumtrapz(Xall(:,i),DX(:,i));
% p4=plot(Xall(:,i),-Q(:,i),'color',colors(i,:));
% , 'LineStyle', 'none','Marker','.');
% hold on

% set(p4,'LineWidth',3)
% set(gca,'FontSize',14)

end
% yline(0)
% xlabel('states (X)')
% 
% ylabel('Potential landscape')

%% potential function, analytical
%  syms r k x l a
% f(x,r,k,l,a)=(r.*x.*(1-x/k)-l*x.^2./(x.^2+a^2));
% Fx=int(f,x);

% xxx=X0:.01:RX0;
% U=(rr*xxx.^2)/2 - lambda*xxx + aa*lambda*atan(xxx/aa) - (rr*xxx.^3)./(3*KK);
% 
% plot(xxx,-U)

%% SK var Md

NL=length(dx)+length(dxML);

for i=1:length(al)
SKL(i)=skewness(Xall(1:NL,i));
SKR(i)=skewness(Xall(NL+1:end,i));

varL(i)=var(Xall(1:NL,i));
varR(i)=var(Xall(NL+1:end,i));

mdL(i)=median(Xall(1:NL,i));
mdR(i)=median(Xall(NL+1:end,i));

end

figure
subplot(3,2,1)

hold on
   for i=1:length(al)
       plot(1-al(i),SKL(i),'Marker','*','LineStyle', 'none')
       
   end
   ylabel('SK left valley')
   xlabel('memory')
subplot(3,2,2)
hold on
   for i=1:length(al)
       plot(1-al(i),SKR(i),'Marker','*','LineStyle', 'none')
   end
      ylabel('SK right valley')
   xlabel('memory')
   
   subplot(3,2,3)

hold on
   for i=1:length(al)
       plot(1-al(i),varL(i),'Marker','*','LineStyle', 'none')
       
   end
   ylabel('var left valley')
   xlabel('memory')
subplot(3,2,4)
hold on
   for i=1:length(al)
       plot(1-al(i),varR(i),'Marker','*','LineStyle', 'none')
   end
      ylabel('VAr right valley')
   xlabel('memory')
   
   subplot(3,2,5)

hold on
   for i=1:length(al)
       plot(1-al(i),mdL(i),'Marker','*','LineStyle', 'none')
       
   end
   ylabel('Median left valley')
   xlabel('memory')
subplot(3,2,6)
hold on
   for i=1:length(al)
       plot(1-al(i),mdR(i),'Marker','*','LineStyle', 'none')
   end
      ylabel('Median right valley')
   xlabel('memory')

%% local min and max
% dP: The prominence of a local minimum (or valley) measures how the valley...
%     stands out with respect to its depth and location relative to other valleys.

for i=1:length(al)
[TF,p] = islocalmin(-Q(:,i));

dP(:,i)=p;
MinU(:,i)=TF;

[TF1,p1] = islocalmax(-Q(:,i));

MaxU(:,i)=TF1;
end


%%
figure

subplot(1,2,1)

hold on
   for i=1:length(al)
       LdP=dP(MinU(:,i),i);
       plot(1-al(i),LdP(1),'Marker','*','LineStyle', 'none')
       
   end
   
subplot(1,2,2)
hold on
   for i=1:length(al)
       LdP=dP(MinU(:,i),i);
       plot(1-al(i),LdP(2),'Marker','*','LineStyle', 'none')
   end
   
%% figure

subplot(1,2,1)

hold on
   for i=1:length(al)
       Lp=-Q(MinU(:,i),i);
       Hp=-Q(MaxU(:,i),i);
       plot([1-al(i),1-al(i)],[Lp(1),Hp(1)])
       
   end
   
subplot(1,2,2)
hold on
   for i=1:length(al)
       Lp=-Q(MinU(:,i),i);
       Hp=-Q(MaxU(:,i),i);
       plot([1-al(i),1-al(i)],[Lp(2),Hp(1)])
   end
   

