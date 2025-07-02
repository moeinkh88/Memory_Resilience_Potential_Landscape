clear
clc

load('PAR.mat')
PAR(:,4)=[];
load('Rdp.mat')
load('CurvR.mat')
% load('p.mat')
load('xmin.mat')
load('xmax.mat')
load('Flat.mat')
load('EngRL.mat')
load('EcoRL.mat')

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
% % % % p(1,II1)=EcoRL(2,II1)-.15;
% % % % 
% % % % for i=1:length(1)
% % % %     for j=1:length(II1)
% % % % tic
% % % % while 1
% % % % 
% % % % [t,x] = FDE_PI2_IM(Order(i),F,JF,t0,T,xmin(II1(j)),h,[PAR(II1(j),:),p(1,II1(j)),10,20]);
% % % % 
% % % % Xmin(j)=x(999);
% % % % 
% % % % if x(end)>xmax(II1(j)) || abs(x(end)-Xmin(j))<1e-5
% % % %     p(1,II1(j))=p(1,II1(j))+0.005;
% % % % elseif x(end)<xmax(II1(j)) || x(end)<1e-5
% % % %     EcoRL(1,II1(j))=(p(1,II1(j))-0.005);    
% % % %     break
% % % % end
% % % % 
% % % % end
% % % % toc    
% % % %  
% % % % 
% % % %     end
% % % % 
% % % % end

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
scatter(Flat(:,1),Rdp(:,1),60,EcoRL(1,:),'filled')
ax = gca;
ax.XDir = 'reverse';
% view(-31,14)
xlabel('Resilience curvature')
ylabel('Resilience depth')

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

% colormap cool
%% plot

figure
scatter(Flat(:,1),EcoRL(1,:),'filled')
% hold on
% scatter(dt, Rdp(:,2),150,EcoRL(2,:),'filled')

ylabel('Ecological resilience')
xlabel('Flatness')

colormap cool
%% plot
figure
scatter(Rdp(:,1),EcoRL(1,:)'./Flat(:,1))
ax = gca;
ax.XDir = 'reverse';
% view(-31,14)
xlabel('Flatness')
ylabel('EngRL/pd')
set(gca,'xscale','log')
set(gca,'yscale','log')
%% plot

figure
scatter(Rdp(:,1),EcoRL(1,:),'filled')

ylabel('Ecological resilience')
xlabel('Potential depth')

ax = gca;
ax.YAxis.Scale ="log";
ax.XAxis.Scale ="log";

colormap cool


%% plot
%colors
Cpr=[0.4940 0.1840 0.5560];
Cgr=[0.4660 0.6740 0.1880];
figure

hold on

indP=find((Rdp(:,11)-Rdp(:,1))>0);
indN=find((Rdp(:,11)-Rdp(:,1))<0);

indP1 = repmat({'Memory increases \newline potential depth'},length(indP),1);
indP2 = repmat({'Memory decreases \newline potential depth'},length(indN),1);

% abs(Rdp(indP,2)-Rdp(indP,1))./Rdp(indP,1);
EcoRdp1=(EcoRL(2,indP)-EcoRL(1,indP))./(EcoRL(2,indP)+EcoRL(1,indP));
% scatter(abs(Rdp(indN,2)-Rdp(indN,1))./Rdp(indN,1),(EcoRL(2,indN)-EcoRL(1,indN))./EcoRL(1,indN),[],Cgr)
EcoRdp2=(EcoRL(2,indN)-EcoRL(1,indN))./(EcoRL(2,indN)+EcoRL(1,indN));

v=violinplot([EcoRdp1';EcoRdp2'],[indP1;indP2],...
'Width',0.5,'Bandwidth',0.02,'ViolinColor',Cgr,'ViolinAlpha',0.004,...
'EdgeColor',[.5,.5,.5],'BoxColor',[.1,.1,.1],'MedianColor',[.1,.1,.1]);

% v.ShowData=true; v.ViolinAlpha=0.1;
 h=v.ScatterPlot; h.MarkerFaceAlpha=1; h.MarkerFaceColor='flat';
 
 h.CData=repmat(h.CData,numel(h.YData),1);
 index=h.YData<3; h.CData(index,1)=Cpr(1); h.CData(index,2)=Cpr(2); h.CData(index,3)=Cpr(3);

ax = gca;
ax.YAxis.Scale ="log";

axis tight

ylabel('Relative effect of memory on ecological resilience')
% xlabel('Memory effects on potential depth')

%% not nice

figure
plot((Rdp(:,11)-Rdp(:,1))./(Rdp(:,11)+Rdp(:,1)),(diff(EcoRL)./sum(EcoRL)),'o')

ylabel('Memory effects on ecological resilience')
xlabel('Memory effects on potential depth')

%% plot

figure
plot(abs(Flat(:,11)-Flat(:,1))./abs(Flat(:,11)+Flat(:,1)),(diff(EcoRL)./sum(EcoRL)),'o')

ylabel('Memory effects on ecological resilience')
xlabel('Memory effects on potential curv')
set(gca,'xscale','log')
set(gca,'yscale','log')


%% 
% g1 = repmat({'EcoRL'},length(EcoRL(1,:)),1);
% 
% figure
% violinplot(diff(EcoRL)./sum(EcoRL),g1);
% ylabel('Impact of memory on EcoRL')

%%
% figure
% scatter(abs(Flat(:,11)-Flat(:,1))./abs(Flat(:,11)+Flat(:,1)),(Rdp(:,11)-Rdp(:,1))./(Rdp(:,11)+Rdp(:,1)),60,(diff(EcoRL)./sum(EcoRL)),'filled')
% ax = gca;
% ax.XDir = 'reverse';
% % view(-31,14)
% xlabel('Memory effects on Flatness')
% ylabel('Memory effects on potential depth')
% 
% cb = colorbar;                                     % create and label the colorbar
% cb.Label.String = 'Relative Effects of Memory on Ecological resilience';
% 
% colormap cool

%% Eco----Curv and pd
figure
p1=scatter(Flat(:,1),Rdp(:,1),60,(EcoRL(1,:)),'filled');
% scatter3(Flat(:,1),Rdp(:,1),(diff(EcoRL)./sum(EcoRL)))
ax = gca;
% ax.XDir = 'reverse';
% view(-31,14)
xlabel('Flatness')
ylabel('Potential depth')
% % 
% ax.YAxis.Scale ="log";
% ax.XAxis.Scale ="log";
% ax.ZAxis.Scale ="log";
set(gca,'ColorScale','log') 
p1.MarkerEdgeColor=[.7,.7,.7];

cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Ecological resilience';

set(gca,'CLim',[0.0520 7.1130]);
colormap(flipud(gray))

% corEcoCurv=corrcoef(Flat(:,1),(EcoRL(1,:)));
% corEcocurv=corEcoCurv(2);
% display(corEcocurv)
% corEcoDp=corrcoef(Rdp(:,1),(EcoRL(1,:)));
% corEcodp=corEcoDp(2);
% display(corEcodp)

%% Eff M Eco----Curv and pd
figure
scatter(Flat(:,1),Rdp(:,1),60,(diff(EcoRL)./sum(EcoRL)),'filled')
% scatter3(Flat(:,1),Rdp(:,1),(diff(EcoRL)./sum(EcoRL)))
% ax = gca;
% ax.XDir = 'reverse';
% view(-31,14)
xlabel('Flatness')
ylabel('Potential depth')
% % 
% ax.YAxis.Scale ="log";
% ax.XAxis.Scale ="log";
% ax.ZAxis.Scale ="log";
% set(gca,'ColorScale','log') 

cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Impact of memory on ecological resilience';

set(gca,'CLim',[0 2.5e-2]);
colormap(flipud(gray))

%% Eff M Eco----Curv

figure
p1=scatter(Flat(:,1),(diff(EcoRL)./sum(EcoRL)),'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7);
% ax = gca;
% ax.XDir = 'reverse';
ylabel('Relative effect of memory on EcoRL')
xlabel('Flatness')
% set(gca,'xscale','log')
set(gca,'yscale','log')
p1.MarkerFaceColor = 'k';
p1.MarkerEdgeColor='w';
% 

corEffMEcoCurv=corrcoef(Flat(:,1),(diff(EcoRL)./sum(EcoRL)));
corEffMEcocurv=corEffMEcoCurv(2);
display(corEffMEcocurv)
%% Eff M Eco---- pd

figure
p11=scatter(Rdp(:,1),(diff(EcoRL)./sum(EcoRL)),'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7);

ylabel('Relative effect of memory on EcoRL')
xlabel('Potential depth')
p11.MarkerFaceColor = 'k';
p11.MarkerEdgeColor='w';

ax = gca;
ax.YAxis.Scale ="log";
ax.XAxis.Scale ="log";

colormap cool

[corEffMEcoDp,p1]=corrcoef(Rdp(:,1),(diff(EcoRL)./sum(EcoRL)));
corEffMEcodp=corEffMEcoDp(2);
display(corEffMEcodp)


%% Eff M Eco---- dt

figure
p11=scatter(dt,(diff(EcoRL)./sum(EcoRL)),'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7);

ylabel('Relative effect of memory on resistance')
xlabel('Distance to basin threshold')
p11.MarkerFaceColor = 'k';
p11.MarkerEdgeColor='w';

ax = gca;
ax.YAxis.Scale ="log";
ax.XAxis.Scale ="log";

colormap cool

[corEffMEcoDp,p1]=corrcoef(dt,(diff(EcoRL)./sum(EcoRL)));
corEffMEcodp=corEffMEcoDp(2);
display(corEffMEcodp)
