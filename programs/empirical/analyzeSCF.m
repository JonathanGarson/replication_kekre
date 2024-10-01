%% analyzeSCF.m
% Plots exposure to corporate profits by wealth quantile from 2004 SCF

tic;

%% Read in dataset

equityshare = xlsread('output/corpexposure_by_wealth.xlsx');

%% Set plotting format
set(0,'defaultaxescolororder',[0 0 0.5]);
set(0,'defaultaxeslinestyleorder',{'-','--','-.',':'});
set(0,'defaultLineLineWidth',2);
set(0,'defaultAxesLineWidth',0.5);
Interpreter = 'latex';

%% Plot

figure;
hax=axes;
plot(equityshare(:,1),equityshare(:,2),'o');
hold on;
axis([1 20 0 0.5]);
hold off;
ax = gca;
ax.XTick = 1:20;
ax.XTickLabel = {'1','','','','','','','','','',...
    '','','','','','','','','','20'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Median of equity share}'],'Interpreter',Interpreter);
xlabel('Quantile of wealth','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat(['output/figures_and_tables/FigureA1']), '-depsc', '-r600','-loose');

disp(['Code completed in ',num2str(toc),' s.']);