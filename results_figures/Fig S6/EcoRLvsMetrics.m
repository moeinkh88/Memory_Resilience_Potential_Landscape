clear
clc


% Define the path to the data folder
dataPath = fullfile('..', '..', 'data');

% Get a list of all .mat files in the data folder
matFiles = dir(fullfile(dataPath, '*.mat'));

% Loop through each .mat file and load it
for i = 1:length(matFiles)
    load(fullfile(dataPath, matFiles(i).name));
end

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



