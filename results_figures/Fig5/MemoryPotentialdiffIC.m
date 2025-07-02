clc
clear

F=@fun1; %logistic model dx/dt=x(1-x)
al=1:-.2:.8; % order of derivatives (memory= 0 and 0.2)
h=0.01; % step size for computing

T=1000; % final time
t0=0; % initial time
X0=.01; % initial condition
RX0=1.6; % initial condition upper the unstable point

for i=1:length(al)
         [t,x] = FDE_PI12_PC(al(i),F,t0,T,X0,h);  % solve equation with and without memory, for the first initial time point=0
    X(i,:)=x;
end
T01=3.62; % second initial time point

indx01=find(abs(t-T01)<10^-6); % find the time close to the second initial time point that there is x value
for i=1:length(al)
         [t1,x1] = FDE_PI12_PC(al(i),F,T01,T,X(i,indx01),h); % solve equation with and without memory, for the second initial time point=T01
    X1(i,:)=x1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% similar evaluations for the upper initial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% value = RX0
for i=1:length(al)
         [Rt,Rx] = FDE_PI12_PC(al(i),F,t0,T,RX0,h); % solve equation with and without memory, for the first initial time point=0
    RX(i,:)=Rx;
end
RT01=.55;  % second initial time point
Rindx01=find(abs(Rt-RT01)<10^-6);
for i=1:length(al)
         [Rt1,Rx1] = FDE_PI12_PC(al(i),F,RT01,T,RX(i,Rindx01),h);% solve equation with and without memory, for the first initial time point=0
    RX1(i,:)=Rx1;
end

%% plotting the dynamics
figure
subplot(1,2,1)
colors=[0 0.4470 0.7410
    0.8500 0.3250 0.0980];

for i=1:length(al)
p1=plot(t,X(i,:),'color',colors(i,:));
hold on
p11=plot(t1,X1(i,:),':');
p11.Color=p1.Color+[.1,.3,0];

p1R=plot(Rt,RX(i,:),'color',colors(i,:));
hold on
p11R=plot(Rt1,RX1(i,:),':');
p11R.Color=p1R.Color+[.1,.3,0];

set(p1,'LineWidth',3)
set(gca,'FontSize',14)
set(p11,'LineWidth',3)
set(gca,'FontSize',14)
set(p1R,'LineWidth',3)
set(gca,'FontSize',14)
set(p11R,'LineWidth',3)
set(gca,'FontSize',14)
end

xxl0=line([T01 T01],[0 X1(1)],'LineStyle',':', 'Color', 'k');
xxl0R=line([RT01 RT01],[0 RX1(1)],'LineStyle',':', 'Color', 'k');

xlabel('Time')
ylabel('States')
axis([0,10,0,RX0])
%% potential landscape

% evaluate the rate of changes based on the derivatives of the solutions
dx=diff(X')./h;
dx1=diff(X1')./h;
dxR=diff(RX')./h;
dx1R=diff(RX1')./h;

subplot(1,2,2)

for i=1:length(al) % evaluating potential landscape
Q1 = cumtrapz(X(i,1:end-1),dx(:,i));
p4=plot(X(i,1:end-1),-Q1,'color',colors(i,:));
hold on
Q11 = cumtrapz(X1(i,1:end-1),dx1(:,i));
p41=plot(X1(i,1:end-1),-Q11-Q1(indx01),':');
p41.Color=p4.Color+[.1,.3,0];

Q1R = cumtrapz(RX(i,1:end-1),dxR(:,i));Q1R=Q1R-abs(Q1R(end)-Q1(end));
p4R=plot(RX(i,1:end-1),-Q1R,'color',colors(i,:));
hold on
Q11R = cumtrapz(RX1(i,1:end-1),dx1R(:,i));
p41R=plot(RX1(i,1:end-1),-Q11R-Q1R(Rindx01),':');
p41R.Color=p4R.Color+[.1,.3,0];

legend('Memory=0','Memory=0.2','With different IC')
set(get(get(p4R,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p41R,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


set(p4,'LineWidth',3)
set(gca,'FontSize',14)
set(p41,'LineWidth',3)
set(gca,'FontSize',14)
set(p4R,'LineWidth',3)
set(gca,'FontSize',14)
set(p41R,'LineWidth',3)
set(gca,'FontSize',14)
end
yyl=yline(0);
set(get(get(yyl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
xxl=line([X1(1) X1(1)], [-.2 -Q1(1)-Q1(indx01)],'LineStyle',':', 'Color', 'k');
set(get(get(xxl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
xxlR=line([RX1(1) RX1(1)], [-.2 -Q11R(1)-Q1R(Rindx01)], 'LineStyle',':','Color', 'k');
set(get(get(xxlR,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% legend('\alpha=1','\alpha=1, X0','\alpha=0.8','\alpha=0.8, X0', '\alpha=0.6','\alpha=0.6, X0')
xlabel('States')
ylabel('Potential energy')
