%% Description
% This script generates the MaxCut plots for Figure 2 in [DYTCU2019].
% You need to run the test file "Test2_MaxCut" first, before this script.
%
% [DYTCU2019] L. Ding, A. Yurtsever, V. Cevher, J.A. Tropp, M. Udell,
% "An Optimal-Storage Approach to Semidefinite Programming using Approximate Complementarity"
% arXiv:1902.03373, 2019.

%%
close all
clearvars

colors = [0 0 0;
    0 0.3 1;
    1 0.1 0.1];

load ./data/G1.mat;
n = size(Problem.A,1);
L = spdiags(Problem.A*ones(n,1),0,n,n) - Problem.A;
C = -0.25*L;
b = ones(n,1);
normb = norm(b,'fro');

load ./results/Test2_MaxCut_Results

optval = res.sol.itr.pobjval;
Xmosek = zeros(n);
[ii,jj] = ndgrid(1:n); % ii and jj are row and column indices respectively
Xmosek(ii>=jj) = res.sol.itr.barx; % fill in values in column-major order
Xmosek = Xmosek + tril(Xmosek,-1).'; % transpose to get result
normXmosek = norm(Xmosek,'fro');

%%
hfig = figure('Position',[100,100,1500,270]);

subplot(141)
myLegend = {};
ind = find(infoAdaGrad_r.y_bar.dualTot - optval <= 0);
loglog(infoAdaGrad_r.iteration(ind),abs(infoAdaGrad_r.y_bar.dualTot(ind) - optval)/abs(optval),'Color',colors(2,:));
hold on
myLegend{end+1} = 'AdaGrad';
ind = find(infoAdaNGD_r.y_bar.dualTot - optval <= 0);
loglog(infoAdaNGD_r.iteration(ind),abs(infoAdaNGD_r.y_bar.dualTot(ind) - optval)/abs(optval),'Color',colors(3,:));
myLegend{end+1} = 'AdaNGD';
ind = find(infoAcceleGrad_r.y_bar.dualTot - optval <= 0);
loglog(infoAcceleGrad_r.iteration(ind),abs(infoAcceleGrad_r.y_bar.dualTot(ind) - optval)/abs(optval),'Color',colors(1,:));
myLegend{end+1} = 'AcceleGrad';
%         loglog(Test1.UniPD100.iteration,abs(Test1.UniPD100.y_bar.dualTot - optval)/abs(optval));
%         myLegend{end+1} = 'UniPD';
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
%         loglog(Test1.UniPD100.iteration,abs(Test1.UniPD100.y_bar.primalObj - optval)/abs(optval));
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
%         loglog(Test1.UniPD100.iteration,Test1.UniPD100.y_bar.primalFeas./normb);
ylabel('feasibility gap','Interpreter','latex','FontSize',18)
grid on; grid minor; grid minor;
subplot(144)
loglog(infoAcceleGrad_r.iteration,infoAcceleGrad_r.y_bar.relDist./normXmosek,'Color',colors(1,:));
hold on
loglog(infoAdaGrad_r.iteration,infoAdaGrad_r.y_bar.relDist./normXmosek,'Color',colors(2,:));
loglog(infoAdaNGD_r.iteration,infoAdaNGD_r.y_bar.relDist./normXmosek,'Color',colors(3,:));
loglog(infoAcceleGrad_3r.iteration,infoAcceleGrad_3r.y_bar.relDist./normXmosek,'--','Color',colors(1,:));
loglog(infoAdaGrad_3r.iteration,infoAdaGrad_3r.y_bar.relDist./normXmosek,'--','Color',colors(2,:));
loglog(infoAdaNGD_3r.iteration,infoAdaNGD_3r.y_bar.relDist./normXmosek,'--','Color',colors(3,:));
%         loglog(Test1.UniPD100.iteration,Test1.UniPD100.y_bar.relDist./normXmosek);hold on
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
    xlim([1,1e5])
end

%% Last edit: Alp Yurtsever - November 6, 2019
