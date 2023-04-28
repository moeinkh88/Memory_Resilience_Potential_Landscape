clear
clc

load('PAR.mat')
PAR(:,4)=[];
load('Rdp.mat')
load('CurvR.mat')
% load('p.mat')
load('xmin.mat')
load('xmax.mat')

%% find Xmin and Xmax

% for i=1:length(PAR)
%     
% xt=0:.0001:11; % steps for polynomial function
% U=-polyint(PAR(i,:));
% U1=polyval(U,xt);
% TF=islocalmin(U1);
% TF1=islocalmax(U1);
% 
% xmin(i)=xt(TF); % critical points: bottom of valley
% xmax(i)=xt(TF1); % critical point: threshold 
% 
% end
%% EcoRL

F=@funPoly;
JF=@JfunPoly;

Order=1:-.2:.8;
t0=0;
T=25;
h=0.01;

% J=1;
% 
% for i=1:length(Order)
%     for j=1:length(PAR)
% tic
%     p1=0.01; p2=10; %initial guess
%     pt0=10;pt1=20;
% while 1
% 
% [t1,x1] = FDE_PI2_IM(Order(i),F,JF,t0,30,xmin(j),h,[PAR(j,:),p1,pt0,pt1]);
% [t2,x2] = FDE_PI2_IM(Order(i),F,JF,t0,30,xmin(j),h,[PAR(j,:),p2,pt0,pt1]);
% if (x1(end)>x1(end-1) && x2(end)>x2(end-1)) || (abs(x1(end)-xmin(j))<1e-4 && abs(x2(end)-xmin(j))<1e-4)
%     p(i,j)=p2;
%     break
%     elseif (x1(end)>x1(end-1) && x2(end)<x2(end-1)) ||  (x1(end)>1e-5 && x2(end)<1e-4)
%     p1=p1*2; p2=p2/2; 
%     if p1>p2
%         p(i,j)=p1;
%     break
%     end
%     elseif (x1(end)<x1(end-1) && x2(end)<x2(end-1)) ||  (x1(end)<1e-5 && x2(end)<1e-4)
% %     else
%     p1=p1/2; p2=p2/2;
% %     if p2<.01 
% % %         pt0=pt0/2;pt1=pt1/2;
% %        I(J)=i; J=1+J; break
% %     end
% end
% 
% end
% toc
%     end
% end
%%
% P(1,:)=p(1,:);
% load('P08.mat')
% P(2,:)=P08;
% % load('EcoRL.mat')
% % p=EcoRL;
% % I1=find(EcoRL>0.1);
% % I2=find(EcoRL<0.1);
% % p(I1)=EcoRL(I1)-0.005;
% % p(I2)=EcoRL(I2)-0.05;
% % 
% % for i=1:length(Order)
% %     for j=1:length(PAR)
% % tic
% % while 1
% % 
% % [t,x] = FDE_PI2_IM(Order(i),F,JF,t0,T,xmin(j),h,[PAR(j,:),p(i,j),10,20]);
% % 
% % Xmin(j)=x(9.99/h);
% % 
% % if x(end)>xmax(j) || abs(x(end)-Xmin(j))<1e-5
% %     p(i,j)=p(i,j)+0.001;
% % elseif x(end)<xmax(j) || x(end)<1e-5
% %     EcoRL(i,j)=(p(i,j)-0.001);    
% %     break
% % end
% % 
% % end
% % toc    
% %  
% % 
% %     end
% % 
% % end

%%
% load('EcoRL.mat')
% I1=find(abs(EcoRL(1,:)-1.243)<.00000001);

% p(1,II)=EcoRL(1,II)/2-1;
p(1,II1)=EcoRL(2,II1)-.15;

for i=1:length(1)
    for j=1:length(II1)
tic
while 1

[t,x] = FDE_PI2_IM(Order(i),F,JF,t0,T,xmin(II1(j)),h,[PAR(II1(j),:),p(1,II1(j)),10,20]);

Xmin(j)=x(999);

if x(end)>xmax(II1(j)) || abs(x(end)-Xmin(j))<1e-5
    p(1,II1(j))=p(1,II1(j))+0.005;
elseif x(end)<xmax(II1(j)) || x(end)<1e-5
    EcoRL(1,II1(j))=(p(1,II1(j))-0.005);    
    break
end

end
toc    
 

    end

end

%%
% for i=1:length(Order)
%     for j=1:length(PAR)
% 
% T=100;
% pt0=10;pt1=20;
% IndxP=pt1/h;
% Eps=1e-3;
% [t,x] = FDE_PI2_IM(Order(i),F,JF,t0,T,xmin(j),h,[PAR(j,:),p(i,j),10,20]);
% while 1
% 
% RX(i,:)=x;
% ind=find(abs(x(IndxP:end)-xmin(j))<Eps);
% 
% if isempty(ind)
%     T=T+500;
%     [t,x] = FDE_PI2_IM(Order(i),F,JF,t0,T,xmin(j),h,[PAR(j,:),p(i,j),10,20]);
% else
% EngRL(i,j)=(2*(abs(xmin(j)-min(x)))/(abs(xmin(j)-min(x))+Eps)-1)/(t(ind(1)));
% break
% end
% end 
% 
%     end
% end

%% plot
dt=xmin-xmax;
figure
% scatter3(CurvR(:,1),dt,Rdp(:,1),40,(EcoRL(2,:)-EcoRL(1,:))./EcoRL(1,:),'filled')
scatter3(CurvR(:,1),dt,Rdp(:,1),100,(EcoRL(2,:)-EcoRL(1,:))./EcoRL(1,:),'filled')
ax = gca;
ax.XDir = 'reverse';
view(-31,14)
xlabel('Curvature')
ylabel('Distance to threshold')
zlabel('Potential depth')

cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Ecological resilience';

colormap cool
%% plot

figure
scatter(dt, Rdp(:,1),150,EcoRL(1,:),'filled')
% hold on
% scatter(dt, Rdp(:,2),150,EcoRL(2,:),'filled')

ylabel('Ecological resilience')
xlabel('Potential depth')

colormap cool
%% plot

figure
scatter(dt, Rdp(:,1),150,EcoRL(2,:)-EcoRL(1,:),'filled')

ylabel('Memory effects on ecological resilience')
xlabel('Potential depth')

colormap cool

%% plot

figure
scatter(Rdp(:,1),EcoRL(2,:)-EcoRL(1,:),'filled')

ylabel('Memory effects on ecological resilience')
xlabel('Potential depth')

colormap cool

%% plot

figure
scatter(Rdp(:,1),EcoRL(1,:),'filled')

ylabel('Ecological resilience')
xlabel('Potential depth')

colormap cool


%% plot

figure
plot((Rdp(:,2)-Rdp(:,1))./Rdp(:,1),(EcoRL(2,:)-EcoRL(1,:))./EcoRL(1,:),'o')

ylabel('Memory effects on ecological resilience')
xlabel('Memory effects on potential depth')


colormap cool

%% plot

figure
plot((CurvR(:,2)-CurvR(:,1))./abs(CurvR(:,1)),(EcoRL(2,:)-EcoRL(1,:))./EcoRL(1,:),'o')

ylabel('Memory effects on ecological resilience')
xlabel('Memory effects on potential curv')


colormap cool
%% plot

figure
% plot(dt,(EcoRL(2,:)-EcoRL(1,:))./EcoRL(1,:),'o')
plot(dt,(EcoRL(2,:)-EcoRL(1,:)),'o')

ylabel('Memory effects on ecological resilience')
xlabel('dt')


colormap cool