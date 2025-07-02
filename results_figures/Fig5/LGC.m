clc
clear
F=@fun1;
al=1;
h=0.01;
% al=.9
    [t,x] = FDE_PI12_PC(al,F,0,10,0.1,h);
    [t1,x1] = FDE_PI12_PC(al,F,0,10,1.5,h);
%%        
colors=[0.4940 0.1840 0.5560];
figure

subplot(1,2,1)
p1=plot(t,x,'color',colors);
set(p1,'LineWidth',3)
set(gca,'FontSize',14)
hold on

p2=plot(t1,x1,'color',colors);
set(p2,'LineWidth',3)
set(gca,'FontSize',14)
xlabel('time')
ylabel('X(t)')

%%potential
T=1.6;
X=0:.01:T;

V= -X.^2./2+X.^3./3;
Xdot=X-X.^2;

subplot(1,2,2)
hold on
p3=plot(X,V,'k');
set(p3,'LineWidth',2)
set(gca,'FontSize',14)
p4=plot(X,Xdot,'--k');
set(p4,'LineWidth',2)
set(gca,'FontSize',14)
ylabel('Potential and derivative function')
xlabel('X(t)')
legend('V','dx/dt','AutoUpdate','off')
yline(0)
xx=1;
yy=linspace(min(V),max(Xdot),40);
plot(xx*ones(size(yy)),yy,':k')
axis([0,T,-.4,max(Xdot)])

h=annotation('arrow', [.73,.78], [.61,.61],'LineWidth',.001);
h1=annotation('arrow', [.9,.8], [.61,.61],'LineWidth',.001);
