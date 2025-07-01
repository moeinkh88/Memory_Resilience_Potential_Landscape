clear
clc

load('Rdp.mat')
load('EngRL.mat')

%%  RelEff=diff(EngRL)/sum(EngRL)
figure
Cpr=[0.4940 0.1840 0.5560];
Cgr=[0.4660 0.6740 0.1880];

hold on

indP=find((Rdp(:,11)-Rdp(:,1))>0);
indN=find((Rdp(:,11)-Rdp(:,1))<0);

indP1 = repmat({'Memory increases \newline potential depth'},length(indP),1);
indP2 = repmat({'Memory decreases \newline potential depth'},length(indN),1);

% abs(Rdp(indP,2)-Rdp(indP,1))./Rdp(indP,1);
EngRdp1=diff(EngRL(:,indP))./sum(EngRL(:,indP));
% scatter(abs(Rdp(indN,2)-Rdp(indN,1))./Rdp(indN,1),(EcoRL(2,indN)-EcoRL(1,indN))./EcoRL(1,indN),[],Cgr)
EngRdp2=diff(EngRL(:,indN))./sum(EngRL(:,indN));

v=violinplot([EngRdp1';EngRdp2'],[indP1;indP2],...
'Width',0.5,'Bandwidth',0.02,'ViolinColor',Cgr,'ViolinAlpha',0.004,...
'EdgeColor',[.5,.5,.5],'BoxColor',[.1,.1,.1],'MedianColor',[.1,.1,.1]);

% v.ShowData=true; v.ViolinAlpha=0.1;
 h=v.ScatterPlot; h.MarkerFaceAlpha=1; h.MarkerFaceColor='flat';
 
 h.CData=repmat(h.CData,numel(h.YData),1);
 index=h.YData<3; h.CData(index,1)=Cpr(1); h.CData(index,2)=Cpr(2); h.CData(index,3)=Cpr(3);

ax = gca;
% ax.YAxis.Scale ="log";

ylabel('Relative impact of memory on engineering resilience')
% xlabel('Memory effects on potential depth')





