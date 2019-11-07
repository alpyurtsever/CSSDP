%% Description
% This script generates the Matrix Completion plots for Figure 2 in [DYTCU2019].
% You need to run the test file "Test2_MatrixCompletion" first, before this script.
%
% [DYTCU2019] L. Ding, A. Yurtsever, V. Cevher, J.A. Tropp, M. Udell,
% "An Optimal-Storage Approach to Semidefinite Programming using Approximate Complementarity"
% arXiv:1902.03373, 2019.

%%

close all
clearvars

load ./data/MatrixCompletionData.mat
load ./results/Test2_MatrixCompletion_Results.mat

normb = norm(Vals);
optval = trace(GTruth);
n = size(GTruth,1);
normGTruth = norm(GTruth,'fro');

%%
hfig = figure('Position',[100,100,1500,270]);
set(hfig,'name','MatrixCompletion','NumberTitle','off');

colors = [0 0 0;
          0 0.3 1;
          1 0.1 0.1];

subplot(141)
myLegend = {};
loglog(infoAdaGrad_3r.iteration,abs(infoAdaGrad_3r.y_bar.dualTot - optval)/abs(optval),'Color',colors(2,:));
hold on
myLegend{end+1} = 'AdaGrad';
loglog(infoAdaNGD_3r.iteration,abs(infoAdaNGD_3r.y_bar.dualTot - optval)/abs(optval),'Color',colors(3,:));
myLegend{end+1} = 'AdaNGD';
loglog(infoAcceleGrad_3r.iteration,abs(infoAcceleGrad_3r.y_bar.dualTot - optval)/abs(optval),'Color',colors(1,:));
myLegend{end+1} = 'AcceleGrad';

hl = legend(myLegend,'Location','SouthWest');
hl.Interpreter = 'latex';
hl.FontSize = 16;

ylabel('dual suboptimality','Interpreter','latex','FontSize',18)
grid on; grid minor; grid minor;
subplot(142)
loglog(infoAcceleGrad_r.iteration,abs(infoAcceleGrad_r.y_bar.primalObj - optval)/abs(optval),'Color',colors(1,:));
hold on
loglog(infoAdaGrad_r.iteration,abs(infoAdaGrad_r.y_bar.primalObj - optval)/abs(optval),'Color',colors(2,:));
loglog(infoAdaNGD_r.iteration,abs(infoAdaNGD_r.y_bar.primalObj - optval)/abs(optval),'Color',colors(3,:));
loglog(infoAcceleGrad_3r.iteration,abs(infoAcceleGrad_3r.y_bar.primalObj - optval)/abs(optval),'--','Color',colors(1,:));
loglog(infoAdaGrad_3r.iteration,abs(infoAdaGrad_3r.y_bar.primalObj - optval)/abs(optval),'--','Color',colors(2,:));
loglog(infoAdaNGD_3r.iteration,abs(infoAdaNGD_3r.y_bar.primalObj - optval)/abs(optval),'--','Color',colors(3,:));
ylabel('primal suboptimality','Interpreter','latex','FontSize',18)
grid on; grid minor; grid minor;
subplot(143)
loglog(infoAcceleGrad_r.iteration,infoAcceleGrad_r.y_bar.primalFeas./normb,'Color',colors(1,:));
hold on
loglog(infoAdaGrad_r.iteration,infoAdaGrad_r.y_bar.primalFeas./normb,'Color',colors(2,:));
loglog(infoAdaNGD_r.iteration,infoAdaNGD_r.y_bar.primalFeas./normb,'Color',colors(3,:));
loglog(infoAcceleGrad_3r.iteration,infoAcceleGrad_3r.y_bar.primalFeas./normb,'--','Color',colors(1,:));
loglog(infoAdaGrad_3r.iteration,infoAdaGrad_3r.y_bar.primalFeas./normb,'--','Color',colors(2,:));
loglog(infoAdaNGD_3r.iteration,infoAdaNGD_3r.y_bar.primalFeas./normb,'--','Color',colors(3,:));
ylabel('feasibility gap','Interpreter','latex','FontSize',18)
grid on; grid minor; grid minor;
subplot(144)
loglog(infoAcceleGrad_r.iteration,infoAcceleGrad_r.y_bar.relDist./normGTruth,'Color',colors(1,:));
hold on
loglog(infoAdaGrad_r.iteration,infoAdaGrad_r.y_bar.relDist./normGTruth,'Color',colors(2,:));
loglog(infoAdaNGD_r.iteration,infoAdaNGD_r.y_bar.relDist./normGTruth,'Color',colors(3,:));
loglog(infoAcceleGrad_3r.iteration,infoAcceleGrad_3r.y_bar.relDist./normGTruth,'--','Color',colors(1,:));
loglog(infoAdaGrad_3r.iteration,infoAdaGrad_3r.y_bar.relDist./normGTruth,'--','Color',colors(2,:));
loglog(infoAdaNGD_3r.iteration,infoAdaNGD_3r.y_bar.relDist./normGTruth,'--','Color',colors(3,:));
ylabel('distance to solution','Interpreter','latex','FontSize',18)
grid on; grid minor; grid minor;


for t = 1:4
    subplot(1,4,t)
    set(gca,'TickLabelInterpreter','latex',...
        'FontSize',15,...
        'XTick',10.^(0:100),...
        'YTick',10.^(-100:100),...
        'TickDir','out',...
        'LineWidth',1,'TickLength',[0.02 0.02]);
    xlabel('iteration','Interpreter','latex','FontSize',18)
    set(findall(gca, 'Type', 'Line' ,'LineStyle', '-'),'LineWidth',2.5);
    set(findall(gca, 'Type', 'Line' ,'LineStyle', '--'),'LineWidth',2.5);
    xlim([1,1e4])
end

%% Last edit: Alp Yurtsever - November 6, 2019
