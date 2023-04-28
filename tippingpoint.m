clc
clear
global rr KK aa h T C

c=.1:.01:.7;
rr=.8;
KK=3;
aa=.2;
t0=0;
h=0.01;
al=1:-.1:1;
X0=-.1;
RX0=3;
F=@funC;
JF=@JfunC;
T=500;

for j=1:length(c)
    C=c(j);
for i=1:length(al)
         [t,x] = FDE_PI2_IM(al(i),F,JF,t0,T,X0,h); 
    X(i,j)=x(end);

         [Rt,Rx] = FDE_PI2_IM(al(i),F,JF,t0,T,RX0,h); 
    RX(i,j)=Rx(end);
end
end

%% plotting

colors=[0.4940 0.1840 0.5560
0.4660 0.6740 0.1880
0.6350 0.0780 0.1840];
figure
for i=1:length(al)
p1=plot(c,X(i,:),'color',colors(i,:));
hold on

p1R=plot(c,RX(i,:),'color',colors(i,:));
legend('Memoryless','Memory=0.1')
set(get(get(p1R,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


set(p1,'LineWidth',3)
set(gca,'FontSize',14)
set(p1R,'LineWidth',3)
set(gca,'FontSize',14)
end
title('Bifurcation diagram')
xlabel('Grazing rate parameter')
ylabel('State')