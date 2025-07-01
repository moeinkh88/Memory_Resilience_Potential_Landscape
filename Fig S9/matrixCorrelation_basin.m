clear
clc

load('Rdp.mat')
load('CurvR.mat')
load('Flat.mat')
load('EngRL.mat')
load('EcoRL.mat')

%% Spearman correlation matrix

% correlation matrix
r = ones(7,7); % Update size to accommodate new factor

% Calculate correlations for the existing factors
EcoEng=corr(EngRL(1,:)',EcoRL(1,:)', 'Type','Spearman');
r(1,2)=EcoEng;

EcoMEco=corr(EcoRL(1,:)',(diff(EcoRL)./sum(EcoRL))','Type','Spearman');
r(1,3)=EcoMEco;

EcoMEng=corr(EcoRL(1,:)',abs(diff(EngRL)./sum(EngRL))','Type','Spearman');%%Absolute
r(1,4)=EcoMEng;

corEcoCurv=corr(Flat(:,1),(EcoRL(1,:))','Type','Spearman');
r(1,5)=corEcoCurv;

corEcoDp=corr(Rdp(:,1),(EcoRL(1,:))','Type','Spearman');
r(1,6)=corEcoDp;

EngMEco=corr(EngRL(1,:)',(diff(EcoRL)./sum(EcoRL))','Type','Spearman');
r(2,3)=EngMEco;

EngMEng=corr(EngRL(1,:)',abs(diff(EngRL)./sum(EngRL))','Type','Spearman');%%Absolute
r(2,4)=EngMEng;

corEngCurv=corr(Flat(:,1),EngRL(1,:)','Type','Spearman');
r(2,5)=corEngCurv;

corEngDp=corr(Rdp(:,1),EngRL(1,:)','Type','Spearman');
r(2,6)=corEngDp;

MEngMEco=corr((diff(EcoRL)./sum(EcoRL))',abs(diff(EngRL)./sum(EngRL))','Type','Spearman');%%Absolute
r(3,4)=MEngMEco;

corEffMEcoCurv=corr(Flat(:,1),(diff(EcoRL)./sum(EcoRL))','Type','Spearman');
r(3,5)=corEffMEcoCurv;

corEffMEcoDp=corr(Rdp(:,1),(diff(EcoRL)./sum(EcoRL))','Type','Spearman');
r(3,6)=corEffMEcoDp;

corEffMEngCurv=corr(Flat(:,1),abs(diff(EngRL)./sum(EngRL))','Type','Spearman');%%Absolute
r(4,5)=corEffMEngCurv;

corEffMEngDp=corr(Rdp(:,1),abs(diff(EngRL)./sum(EngRL))','Type','Spearman');%%Absolute
r(4,6)=corEffMEngDp;

dpCurv=corr(Rdp(:,1),Flat(:,1),'Type','Spearman');
r(5,6)=dpCurv;

% Normalize flatness and depth
normalized_flatness = (abs(CurvR(:,1)) - min(abs(CurvR(:,1)))) / (max(abs(CurvR(:,1))) - min(abs(CurvR(:,1))));
normalized_depth = (Rdp(:,1) - min(Rdp(:,1))) / (max(Rdp(:,1)) - min(Rdp(:,1)));

% Calculate the new metric M (Basin Index)
M = 1/2 * normalized_flatness + 1/2 * normalized_depth;

% Calculate correlations involving M (Basin Index)
corMEcoRL = corr(M, EcoRL(1,:)', 'Type', 'Spearman');
r(1,7) = corMEcoRL;
r(7,1) = corMEcoRL;

corMEngRL = corr(M, EngRL(1,:)', 'Type', 'Spearman');
r(2,7) = corMEngRL;
r(7,2) = corMEngRL;

corMEffMEco = corr(M, (diff(EcoRL)./sum(EcoRL))', 'Type', 'Spearman');
r(3,7) = corMEffMEco;
r(7,3) = corMEffMEco;

corMEffMEng = corr(M, abs(diff(EngRL)./sum(EngRL))', 'Type', 'Spearman');
%%Not Absolute
r(4,7) = corMEffMEng;
r(7,4) = corMEffMEng;

corMFlat = corr(M, Flat(:,1), 'Type', 'Spearman');
r(5,7) = corMFlat;
r(7,5) = corMFlat;

corMDp = corr(M, Rdp(:,1), 'Type', 'Spearman');
r(6,7) = corMDp;
r(7,6) = corMDp;

% Update figure
figure
% labels
labels={'Resistance','Resilience','Relative impact of \newline Memory on Resistance','Relative impact of \newline Memory on Resilience','Flatness','Depth', 'Basin Index'};
% scatter plot
n = size(r, 1);
y = triu(repmat(n+1, n, n) - (1:n)') + 0.5;
x = triu(repmat(1:n, n, 1)) + 0.5;
x(x == 0.5) = NaN;
scatter(x(:), y(:), 1500.*abs(r(:)), r(:), 's', 'filled')

% Convert numeric labels to strings with two decimal places
labelStrings = arrayfun(@(x) sprintf('%.2f', x), r, 'UniformOutput', false);

for i = 2:length(x(:))
    if mod(i, 8) == 1
       text(x(i), y(i), {}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k');
    else
        text(x(i), y(i), labelStrings{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k');
    end
end

% Show labels
text(1:n, (n:-1:1) + 0.5, labels, 'HorizontalAlignment', 'right')
text((1:n) + 0.5, repmat(n + 1, n, 1), labels, 'HorizontalAlignment', 'left', 'Rotation', -270)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nColors = 256;
softness = 0.3;  % lower = more pastel (0.5~0.8 looks good)

% Pastel blue to white
blue = [linspace(softness,1,nColors/2)', linspace(softness,1,nColors/2)', ones(nColors/2,1)];

% White to pastel red
red  = [ones(nColors/2,1), linspace(1,softness,nColors/2)', linspace(1,softness,nColors/2)'];

cmap = [blue; red];
colormap(gca, cmap);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = gca;
colorbar(h);
h.Visible = 'off';
h.Position(4) = h.Position(4)*0.8;
axis(h, 'equal')
caxis([-1, 1]);
cb = colorbar(h);

annotation('textbox', cb.Position, 'FitBoxToText', 'off', 'EdgeColor', 'none');

title('(b) Spearman correlation')