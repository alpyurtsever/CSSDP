%% Description
% This script generates the plots for Figure 1 in [DYTCU2019].
% You need to run the test file "Test1_MaxCut" first, before this script.
%
% [DYTCU2019] L. Ding, A. Yurtsever, V. Cevher, J.A. Tropp, M. Udell,
% "An Optimal-Storage Approach to Semidefinite Programming using Approximate Complementarity"
% arXiv:1902.03373, 2019.

%%
clearvars
close all

load ./data/G1.mat;
n = size(Problem.A,1);
L = spdiags(Problem.A*ones(n,1),0,n,n) - Problem.A;
C = -0.25*L;
b = ones(n,1);

load ./results/Test1_MaxCut_Results.mat

optval = res.sol.itr.pobjval;
Xmosek = zeros(n);
[ii,jj] = ndgrid(1:n); % ii and jj are row and column indices respectively
Xmosek(ii>=jj) = res.sol.itr.barx; % fill in values in column-major order
Xmosek = Xmosek + tril(Xmosek,-1).'; % transpose to get result
%%% Dual norm from Mosek for rough tuning!
yDual = res.sol.itr.y;

%% MinFeas
hfig{1} = figure('Position',[100,100,1100,270]);
set(hfig{1},'name','MaxCutMinFeasTest','NumberTitle','off');
alpha = 800;


optDual = yDual'*b;
relDualError = abs(Out.dualObj+alpha*Out.dualFeas - optDual)/abs(optDual);
relDualError_mean = mean(relDualError,2);
relDualError_min = min(relDualError,[],2);
relDualError_max = max(relDualError,[],2);
MF_PrimalSuboptimality_r = abs(Out.MinFeas.r1.primalObj-optval)/abs(optval);
MF_PrimalSuboptimality_3r = abs(Out.MinFeas.r3.primalObj-optval)/abs(optval);
MF_FeasibiltyGap_r = Out.MinFeas.r1.primalFeas/norm(b);
MF_FeasibiltyGap_3r = Out.MinFeas.r3.primalFeas/norm(b);
MF_Dist2Solution_r = Out.MinFeas.r1.relDist/norm(Xmosek,'fro');
MF_Dist2Solution_3r = Out.MinFeas.r3.relDist/norm(Xmosek,'fro');
MO_PrimalSuboptimality_r = abs(Out.MinObj.r1.primalObj-optval)/abs(optval);
MO_PrimalSuboptimality_3r = abs(Out.MinObj.r3.primalObj-optval)/abs(optval);
MO_FeasibiltyGap_r = Out.MinObj.r1.primalFeas/norm(b);
MO_FeasibiltyGap_3r = Out.MinObj.r3.primalFeas/norm(b);
MO_Dist2Solution_r = Out.MinObj.r1.relDist/norm(Xmosek,'fro');
MO_Dist2Solution_3r = Out.MinObj.r3.relDist/norm(Xmosek,'fro');
subplot(131)

hold on
scatter(relDualError(:),MF_PrimalSuboptimality_r(:),'b','Marker','x');
scatter(relDualError(:),MF_PrimalSuboptimality_3r(:),'r','Marker','o');
plot(relDualError_mean,mean(abs(Out.MinFeas.r1.primalObj-optval)/abs(optval),2),'b')
plot(relDualError_mean,mean(abs(Out.MinFeas.r3.primalObj-optval)/abs(optval),2),'r--')
ylim([0.9999e-7,2])
ylabel('primal suboptimality','Interpreter','latex','FontSize',18)
xlim([1e-7,1])
xlabel('dual suboptimality','Interpreter','latex','FontSize',18)
title('MinFeasSDP','Interpreter','latex','FontSize',14)

subplot(132)
hold on
scatter(relDualError(:),MF_FeasibiltyGap_r(:),'b','Marker','x');
scatter(relDualError(:),MF_FeasibiltyGap_3r(:),'r','Marker','o');
plot(relDualError_mean,mean(Out.MinFeas.r1.primalFeas,2)/norm(b),'b')
plot(relDualError_mean,mean(Out.MinFeas.r3.primalFeas,2)/norm(b),'r--')
ylim([0.9999e-5,2])
ylabel('feasibility gap','Interpreter','latex','FontSize',18)
xlim([1e-7,1])
xlabel('dual suboptimality','Interpreter','latex','FontSize',18)
title('MinFeasSDP','Interpreter','latex','FontSize',14)

subplot(133)
hold on
scatter(relDualError(:),MF_Dist2Solution_r(:),'b','Marker','x');
scatter(relDualError(:),MF_Dist2Solution_3r(:),'r','Marker','o');
plot(relDualError_mean,mean(Out.MinFeas.r1.relDist,2)/norm(Xmosek,'fro'),'b')
plot(relDualError_mean,mean(Out.MinFeas.r3.relDist,2)/norm(Xmosek,'fro'),'r--')
ylim([0.9999e-4,2])
ylabel('distance to solution','Interpreter','latex','FontSize',18)
xlim([1e-7,1])
xlabel('dual suboptimality','Interpreter','latex','FontSize',18)
title('MinFeasSDP','Interpreter','latex','FontSize',14)

for t = 1:3
    subplot(1,3,t)
    set(gca,'TickLabelInterpreter','latex',...
        'FontSize',15,...
        'XTick',10.^(-100:100),...
        'YTick',10.^(-100:100),...
        'TickDir','out',...
        'LineWidth',1,'TickLength',[0.02 0.02]);
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    set(gca, 'XDir','reverse', 'XScale', 'log', 'YScale','log')
    clearvars alpha
    alpha(0.15);
    grid on; grid minor; grid minor;
    box on
end

%% MinObj
hfig{2} = figure('Position',[100,100,1100,270]);
set(hfig{2},'name','MaxCutMinObjTest','NumberTitle','off');

subplot(131)
hold on
scatter(relDualError(:),MO_PrimalSuboptimality_r(:),'b','Marker','x');
scatter(relDualError(:),MO_PrimalSuboptimality_3r(:),'r','Marker','o');
plot(relDualError_mean,mean(abs(Out.MinObj.r1.primalObj-optval)/abs(optval),2),'b')
plot(relDualError_mean,mean(abs(Out.MinObj.r3.primalObj-optval)/abs(optval),2),'r--')
ylim([0.9999e-7,2])
ylabel('primal suboptimality','Interpreter','latex','FontSize',18)
xlim([1e-7,1])
xlabel('dual suboptimality','Interpreter','latex','FontSize',18)
title('MinObjSDP','Interpreter','latex','FontSize',14)

subplot(132)
hold on
scatter(relDualError(:),MO_FeasibiltyGap_r(:),'b','Marker','x');
scatter(relDualError(:),MO_FeasibiltyGap_3r(:),'r','Marker','o');
plot(relDualError_mean,mean(Out.MinObj.r1.primalFeas,2)/norm(b),'b')
plot(relDualError_mean,mean(Out.MinObj.r3.primalFeas,2)/norm(b),'r--')
ylim([0.9999e-5,2])
ylabel('feasibility gap','Interpreter','latex','FontSize',18)
xlim([1e-7,1])
xlabel('dual suboptimality','Interpreter','latex','FontSize',18)
title('MinObjSDP','Interpreter','latex','FontSize',14)

subplot(133)
hold on
plot(relDualError_mean,mean(Out.MinObj.r1.relDist,2)/norm(Xmosek,'fro'),'b')
plot(relDualError_mean,mean(Out.MinObj.r3.relDist,2)/norm(Xmosek,'fro'),'r--')
scatter(relDualError(:),MO_Dist2Solution_r(:),'b','Marker','x');
scatter(relDualError(:),MO_Dist2Solution_3r(:),'r','Marker','o');
ylim([0.9999e-4,2])
ylabel('distance to solution','Interpreter','latex','FontSize',18)
xlim([1e-7,1])
xlabel('dual suboptimality','Interpreter','latex','FontSize',18)
title('MinObjSDP','Interpreter','latex','FontSize',14)
hl = legend('~$r=r^\star$','~$r = 3r^\star$');
hl.Interpreter = 'latex';
hl.FontSize = 17;
hl.Location = 'SouthWest';

for t = 1:3
    subplot(1,3,t)
    set(gca,'TickLabelInterpreter','latex',...
        'FontSize',15,...
        'XTick',10.^(-100:100),...
        'YTick',10.^(-100:100),...
        'TickDir','out',...
        'LineWidth',1,'TickLength',[0.02 0.02]);
    set(findall(gca, 'Type', 'Line' ),'LineWidth',2);
    set(gca, 'XDir','reverse', 'XScale', 'log', 'YScale','log')
    clearvars alpha
    alpha(0.15);
    grid on; grid minor; grid minor;
    box on
end

%% Last edit: Alp Yurtsever - November 6, 2019
