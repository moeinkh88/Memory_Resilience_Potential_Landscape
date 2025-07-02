clc
clear


% Define the path to the data folder
dataPath = fullfile('..', '..', 'data');

% Get a list of all .mat files in the data folder
matFiles = dir(fullfile(dataPath, '*.mat'));

% Loop through each .mat file and load it
for i = 1:length(matFiles)
    load(fullfile(dataPath, matFiles(i).name));
end

xt=0:.00001:5; % steps for polynomial function
Order=1:-.02:.8; % order of derivatives (1-memory) 

Nsample=1000;
for i=1:Nsample

%FLatness
% Flat(i,:)=ceil(max(max(abs(CurvR))))+CurvR(i,:);

[R,pv2]=corrcoef(Rdp(i,:),1-Order); % correlation between Depth and memory
corrRdp(i)=R(2);
[R,pv7]=corrcoef(Flat(i,:),1-Order); % correlation between Flatness and memory
corrCurvR(i)=R(2);

[RHO,pv9] = corr(Flat(i,:)',1-Order','Type','Spearman'); % correlation between Depth and memory (Spearman)
corrCurvRSpear(i)=RHO;


pv(i)=mean([pv2(2),pv7(2),pv9]); %average of p values

end

%colors
Cpr=[0.4940 0.1840 0.5560];
Cgr=[0.4660 0.6740 0.1880];

indPc=find(corrRdp>0); % memory increases the depth
indNc=find(corrRdp<0); % memory decreases the depth

%% Fig a: correlation depth and memory
% note: the dots inside the violin plots are randomly chosen

figure


indP1 = repmat({'Memory increases \newline potential depth'},length(indPc),1);
indP2 = repmat({'Memory decreases \newline potential depth'},length(indNc),1);

% v=violinplot(corrRdp,[], ...
v=violinplot([corrRdp(indPc)';corrRdp(indNc)'],[indP1;indP2],...
'Width',0.5,'Bandwidth',0.04,'ViolinColor',Cgr,'ViolinAlpha',0.004,...
'EdgeColor',[.5,.5,.5],'BoxColor',[.1,.1,.1],'MedianColor',[.1,.1,.1]);

% v.ShowData=true; v.ViolinAlpha=0.1;
 h=v.ScatterPlot; h.MarkerFaceAlpha=1; h.MarkerFaceColor='flat';
 
 h.CData=repmat(h.CData,numel(h.YData),1);
 index=h.YData<0; h.CData(index,1)=Cpr(1); h.CData(index,2)=Cpr(2); h.CData(index,3)=Cpr(3);
ylabel('Correlation: memory and potential depth')
% set(gca,'FontSize',14)

%% Fig b: Effects of memory on curvature
% note: the dots inside the violin plots are randomly chosen

figure

indP11 = repmat({'Pearson'},length(corrCurvR),1);
indP22 = repmat({'Spearman'},length(corrCurvRSpear),1);


violinplot([corrCurvR;corrCurvRSpear],[indP11';indP22'],...
'Width',0.5,'Bandwidth',0.01,'ViolinColor',[.5,.5,.5],'ViolinAlpha',0.004,...
'GroupOrder',{'Spearman','Pearson'}, ...
'EdgeColor',[.5,.5,.5],'BoxColor',[.1,.1,.1],'MedianColor',[.1,.1,.1]);
ylabel('Correlation: memory and flatness')
% set(gca,'FontSize',14)

ylim([0.5 1])

%% Fig c-d: dynamics of some samples


figure
pl1(1)=plot(0:.001:2.265, Potential{indPc(5)}(:,1),'LineWidth',4,'Color','k');
hold on
for i = 2:11
    pl1(i)=plot(0:.001:2.265, Potential{indPc(5)}(:,i), 'Color', Cgr + 0.055*(i-2)-.18,'LineWidth',4);
end
axis tight

strings={','};
% Generate legend or title for showing memory
k =1; ind = 1; 
while k < length(Order)+1
   set(pl1(ind),'DisplayName', [num2str(1-Order(k), '%.2f') ' ' strings{1}, ...
       num2str(Rdp(indPc(5),k), '%.2f') ' ' strings{1}, num2str(Flat(indPc(5),k), '%.2f')])
   k = k+1; ind = ind+1;
end

hleg1=legend('show');
title(hleg1,'Memory, Dept, Flat:')
set(hleg1,'Location','bestoutside','FontSize',10)
xlabel('States')
ylabel('Potential energy')

figure

pl2(1)=plot(0:.002:2.565, Potential{indNc(25)}(:,1),'LineWidth',3,'Color','k');
hold on
for i = 2:11
    pl2(i)=plot(0:.002:2.565, Potential{indNc(25)}(:,i), 'Color', Cpr + 0.067*(i-2)-.18,'LineWidth',3);
end
axis tight

% Generate legend or title for showing memory
k =1; ind = 1; 
while k < length(Order)+1
      set(pl2(ind),'DisplayName', [num2str(1-Order(k), '%.2f') ' ' strings{1}, ...
       num2str(Rdp(indNc(25),k), '%.2f') ' ' strings{1}, num2str(Flat(indNc(25),k), '%.2f')])
   k = k+1; ind = ind+1;
end
hleg1=legend('show');
title(hleg1,'Memory, Dept, Flat:')
set(hleg1,'Location','bestoutside','FontSize',10)
xlabel('States')
ylabel('Potential energy')


