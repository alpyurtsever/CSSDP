%% Description
% This script implements the MaxCut experiment to show that our approach is
% "Compatible with many dual solvers." See Figure 2 in [DYTCU2019].
% This experiment uses the "G1" dataset from the link:
% "https://www.cise.ufl.edu/research/sparse/mat/Gset/G1.mat". 
% The first time you run the code, it will download the dataset under
% the "./data" folder. Please read the README file for the dataset at:
% "https://www.cise.ufl.edu/research/sparse/mat/Gset/README.txt"
% By running this script, you agree to download the dataset and its terms
% and conditions of use. Run "PlotFig2_MaxCut" to generate the plots
% once the experiments are finished.
%
% [DYTCU2019] L. Ding, A. Yurtsever, V. Cevher, J.A. Tropp, M. Udell,
% "An Optimal-Storage Approach to Semidefinite Programming using Approximate Complementarity"
% arXiv:1902.03373, 2019.

%% Beginning of the code
clearvars;
if ~exist('./data','dir'), mkdir data; end
if ~exist('./data/G1.mat','file')
    websave('./data/G1.mat','https://www.cise.ufl.edu/research/sparse/mat/Gset/G1.mat');
    clearvars;
end

%%% fix the random seed
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

clearvars Problem L;

optval = res.sol.itr.pobjval;
Xmosek = zeros(n);
[ii,jj] = ndgrid(1:n); % ii and jj are row and column indices respectively
Xmosek(ii>=jj) = res.sol.itr.barx; % fill in values in column-major order
Xmosek = Xmosek + tril(Xmosek,-1).'; % transpose to get result

%%% setup parameter
singularValues = sort(max(eig(Xmosek),0),'descend');
r = sum(singularValues >= singularValues(1)./1e3);
fprintf('Approximate solution rank is %d.\n',r);

%%% solve the dual problem
alpha = 1.1*n;
Aop = @(V,S) sum((V*S).*V,2);
Atop = @(y,V) repmat(y,1,size(V,2)).*V;
err{1} = 'relDist';
err{2} = @(y,S,V) norm(Xmosek-V*S*V','fro');

maxit = 1e5; % number of iteration

%% NOTE about choosing 'D' below
% We tuned all methods by trying the powers of 10 for 'D', and set here the
% best performing value.

%% AcceleGrad
[~,infoAcceleGrad_r] = AcceleGradSDP( C, b, Aop, Atop, alpha, r, 'D', 100, 'errfncs', err, 'maxit', maxit);

[~,infoAcceleGrad_3r] = AcceleGradSDP( C, b, Aop, Atop, alpha, 3*r, 'D', 100, 'errfncs', err, 'maxit', maxit);

%% AdaNGD
[~,infoAdaNGD_r] = AdaNGDSDP( C, b, Aop, Atop, alpha, r, 1.1, 'D', 100, 'errfncs', err, 'maxit', maxit);

[~,infoAdaNGD_3r] = AdaNGDSDP( C, b, Aop, Atop, alpha, 3*r, 1.1, 'D', 100, 'errfncs', err, 'maxit', maxit);

%% AdaGrad
[~,infoAdaGrad_r] = AdaGradSDP( C, b, Aop, Atop, alpha, r, 'D', 10, 'errfncs', err, 'maxit', maxit);

[~,infoAdaGrad_3r] = AdaGradSDP( C, b, Aop, Atop, alpha, 3*r, 'D', 10, 'errfncs', err, 'maxit', maxit);

%% Save Results
if ~exist('./results','dir'), mkdir results; end
save('./results/Test2_MaxCut_Results.mat','info*','res');

%% Last edit: Alp Yurtsever - November 6, 2019
