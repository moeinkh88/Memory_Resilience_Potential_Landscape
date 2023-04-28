clear
clc

load('PAR.mat')
PAR(:,4)=[];
load('Rdp.mat')
load('CurvR.mat')
% load('p.mat')
load('xmin.mat')
load('xmax.mat')
load('EcoRL.mat')
%% EcoRL

F=@funPoly;
JF=@JfunPoly;

Order=1:-.2:.8;
t0=0;
T=100;
h=0.01;

p=EcoRL./10;
%%
for i=1:length(Order)
    for j=1:length(PAR)
tic
T=100;
pt0=10;pt1=20;
IndxP=pt1/h;
Eps=1e-3;
[t,x] = FDE_PI2_IM(Order(i),F,JF,t0,T,xmin(j),h,[PAR(j,:),p(i,j),10,20]);
while 1

ind=find(abs(x(IndxP:end)-xmin(j))<Eps);

if isempty(ind)
    T=T+500;
    [t,x] = FDE_PI2_IM(Order(i),F,JF,t0,T,xmin(j),h,[PAR(j,:),p(i,j),10,20]);
else
EngRL(i,j)=(2*(abs(xmin(j)-min(x)))/(abs(xmin(j)-min(x))+Eps)-1)/(t(ind(1)));
break
end
end 
toc
    end
end

%% plot

dt=xmin-xmax;
figure
% scatter3(CurvR(:,1),dt,Rdp(:,1),40,(EcoRL(2,:)-EcoRL(1,:))./EcoRL(1,:),'filled')
scatter3(CurvR(:,1),dt,Rdp(:,1),60,(EngRL(1,:)-EngRL(2,:)),'filled')
ax = gca;
ax.XDir = 'reverse';
view(-31,14)
xlabel('Curvature')
ylabel('Distance to threshold')
zlabel('Potential depth')

cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Engineering resilience';

colormap cool