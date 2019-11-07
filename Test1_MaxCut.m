%% Description
% This script implements the synthetic MaxCut experiment. See Figure 1 in
% [DYTCU2019]. This experiment uses the "G1" dataset from the link:
% "https://www.cise.ufl.edu/research/sparse/mat/Gset/G1.mat". 
% The first time you run the code, it will download the dataset under
% the "./data" folder. Please read the README file for the dataset at:
% "https://www.cise.ufl.edu/research/sparse/mat/Gset/README.txt"
% By running this script, you agree to download the dataset and its terms
% and conditions of use. Run "PlotFig1_MaxCut" to generate the plots
% once the experiments are finished.
%
% [DYTCU2019] L. Ding, A. Yurtsever, V. Cevher, J.A. Tropp, M. Udell,
% "An Optimal-Storage Approach to Semidefinite Programming using Approximate Complementarity"
% arXiv:1902.03373, 2019.

%% Beginning of the code
clearvars;

addpath methods

if ~exist('./data','dir'), mkdir data; end
if ~exist('./data/G1.mat','file')
    websave('./data/G1.mat','https://www.cise.ufl.edu/research/sparse/mat/Gset/G1.mat');
end

rng(0,'twister');

load data/G1.mat;
n = size(Problem.A,1);
L = spdiags(Problem.A*ones(n,1),0,n,n) - Problem.A;
C = -0.25*L;
b = ones(n,1);

% Solve with MoSeK for ground truth
prob.bardim = size(L,1);
Ltril = tril(L);
clearvars L;
[prob.barc.subk, prob.barc.subl, prob.barc.val] = find(-0.25.*Ltril);
prob.barc.subj  = ones(size(prob.barc.subk));
prob.a          = sparse(n,n);
prob.bara.subi 	= (1:n)';
prob.bara.subj 	= ones(n,1);
prob.bara.subk 	= (1:n)';
prob.bara.subl 	= (1:n)';
prob.bara.val  	= ones(size(prob.bara.subi));
prob.blc 		= ones(size(prob.bara.subi));
prob.buc 		= ones(size(prob.bara.subi));
[~,res]         = mosekopt('minimize',prob);

optval = res.sol.itr.pobjval;
Xmosek = zeros(n);
[ii,jj] = ndgrid(1:n); % ii and jj are row and column indices respectively
Xmosek(ii>=jj) = res.sol.itr.barx; % fill in values in column-major order
Xmosek = Xmosek + tril(Xmosek,-1).'; % transpose to get result
%%% Dual norm from Mosek for rough tuning!
yDual = res.sol.itr.y;

%%% setup parameter 
singularValues = sort(max(eig(Xmosek),0),'descend');
r1 = sum(singularValues >= singularValues(1)./1e3);
r3 = 3*r1;
info.r1 = r1;
info.r3 = r3;

%%% solve the dual problem 
alpha = 1.1*n;
info.alpha = alpha;
Aop = @(V,S) sum((V*S).*V,2);
Atop = @(y,V) repmat(y,1,size(V,2)).*V;

%%% Recover primal

MCiter = 40;
noiseSweep = [0,10.^(-6.5:0.5:1)];
Out.dualObj = nan(length(noiseSweep),MCiter);
Out.dualFeas = nan(length(noiseSweep),MCiter);
Out.MinFeas.r1.primalObj = nan(length(noiseSweep),MCiter);
Out.MinFeas.r1.primalFeas = nan(length(noiseSweep),MCiter);
Out.MinFeas.r1.relDist = nan(length(noiseSweep),MCiter);
Out.MinFeas.r3.primalObj = nan(length(noiseSweep),MCiter);
Out.MinFeas.r3.primalFeas = nan(length(noiseSweep),MCiter);
Out.MinFeas.r3.relDist = nan(length(noiseSweep),MCiter);
Out.MinObj.r1.primalObj = nan(length(noiseSweep),MCiter);
Out.MinObj.r1.primalFeas = nan(length(noiseSweep),MCiter);
Out.MinObj.r1.relDist = nan(length(noiseSweep),MCiter);
Out.MinObj.r3.primalObj = nan(length(noiseSweep),MCiter);
Out.MinObj.r3.primalFeas = nan(length(noiseSweep),MCiter);
Out.MinObj.r3.relDist = nan(length(noiseSweep),MCiter);

optsEigs.p = 2*r3;
optsEigs.isreal = 1;
optsEigs.issym = issymmetric(C);

for nl = 1:length(noiseSweep)
    nn = noiseSweep(nl);
    for MC = 1:MCiter
        yy = yDual + nn*normc(randn(size(yDual)))*norm(yDual);
        eigsOp = @(x) C*x - Atop(yy,x);
        
        [Vr3,Dr3] = eigs(eigsOp ,n ,r3 ,'sa', optsEigs);
        [~,ind] = sort(diag(Dr3),'ascend');
        Vr1 = Vr3(:,ind(1:r1));
        Dr1 = Dr3(ind(1:r1),ind(1:r1));
        %%% for r = rstar
        
        cvx_begin
        cvx_quiet true
        variable Sr1(r1,r1) semidefinite
        minimize norm(Aop(Vr1,Sr1)-b,'fro')
        cvx_end
        delta = cvx_optval;

        fprintf('r1: noiseLevel %d/%d - MonteCarlo %d \n', ...
            nl, length(noiseSweep), MC);

        Out.MinFeas.r1.primalObj(nl,MC) = sum(sum((C*Vr1).*(Vr1*Sr1)));
        Out.MinFeas.r1.primalFeas(nl,MC) = norm(Aop(Vr1,Sr1) - b,'fro');
        Out.MinFeas.r1.relDist(nl,MC) = norm(Xmosek-Vr1*Sr1*Vr1','fro');

        cvx_begin
        cvx_quiet true
        variable Sr1(r1,r1) semidefinite
        minimize trace((Vr1'*C*Vr1)*Sr1)
        subject to;
        norm(Aop(Vr1,Sr1) - b) <= 1.1*delta;
        cvx_end
        
        Out.MinObj.r1.primalObj(nl,MC) = sum(sum((C*Vr1).*(Vr1*Sr1)));
        Out.MinObj.r1.primalFeas(nl,MC) = norm(Aop(Vr1,Sr1) - b,'fro');
        Out.MinObj.r1.relDist(nl,MC) = norm(Xmosek-Vr1*Sr1*Vr1','fro');
        
        %%% for r = 3rstar
        cvx_begin
        cvx_quiet true
        variable Sr3(r3,r3) semidefinite
        minimize  norm(Aop(Vr3,Sr3)-b,'fro')
        cvx_end
        delta = cvx_optval;
        
        fprintf('r3: noiseLevel %d/%d - MonteCarlo %d \n', ...
            nl, length(noiseSweep), MC);
        Out.MinFeas.r3.primalObj(nl,MC) = sum(sum((C*Vr3).*(Vr3*Sr3)));
        Out.MinFeas.r3.primalFeas(nl,MC) = norm(Aop(Vr3,Sr3) - b,'fro');
        Out.MinFeas.r3.relDist(nl,MC) = norm(Xmosek-Vr3*Sr3*Vr3','fro');
        
        cvx_begin
        cvx_quiet true
        variable Sr3(r3,r3) semidefinite
        minimize trace((Vr3'*C*Vr3)*Sr3)
        subject to;
        norm(Aop(Vr3,Sr3) - b) <= 1.1*delta;
        cvx_end
        
        Out.MinObj.r3.primalObj(nl,MC) = sum(sum((C*Vr3).*(Vr3*Sr3)));
        Out.MinObj.r3.primalFeas(nl,MC) = norm(Aop(Vr3,Sr3) - b,'fro');
        Out.MinObj.r3.relDist(nl,MC) = norm(Xmosek-Vr3*Sr3*Vr3','fro');
        
        Out.dualObj(nl,MC) = yy'*b;
        Out.dualFeas(nl,MC) = abs(min(min(diag(Dr3)),0));
    end
    
end

if ~exist('./results','dir'), mkdir results; end
save('./results/Test1_MaxCut_Results.mat','Out','res');

%% Last edit: Alp Yurtsever - November 6, 2019
