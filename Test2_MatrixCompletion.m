%% Description
% This script implements the Matrix Completion experiment to show that our 
% approach is "Compatible with many dual solvers." 
% Run "PlotFig2_MatrixCompletion" to generate the plots once the 
% experiments are finished.
%
% [DYTCU2019] L. Ding, A. Yurtsever, V. Cevher, J.A. Tropp, M. Udell,
% "An Optimal-Storage Approach to Semidefinite Programming using Approximate Complementarity"
% arXiv:1902.03373, 2019.

%% Beginning of the code
clearvars;
addpath methods
addpath utils

%%% fix the random seed
rng(0,'twister'); 

load('./data/MatrixCompletionData.mat')

n = size(GTruth,1);
m = length(Vals);
C = speye(n);
b = Vals;

%%% solve the dual problem
Aop = @(V,Z) sum( (V(SubI,:)*Z).*V(SubJ,:),2);  % generalized version, a bit more efficient when param.rank is high!
% Atop = @(y,x) sparse(SubI,SubJ,y,n,n,m)*x;
Atop = @(y,x) spmult(SubI,SubJ,y,n,x); % Faster then generating a sparse matrix again and again!

err = {};
err{1} = 'relDist';
err{2} = @(y,S,V) norm(GTruth - V*S*V','fro');

maxit = 1e4; % number of iteration

%% NOTE about choosing 'D' below
% We tuned all methods by trying the powers of 10 for 'D', and set here the
% best performing value.

%% AcceleGrad
[~,infoAcceleGrad_r] = AcceleGradSDP( C, b, Aop, Atop, alpha, r, 'D', 10, 'errfncs', err, 'maxit', maxit);

[~,infoAcceleGrad_3r] = AcceleGradSDP( C, b, Aop, Atop, alpha, 3*r, 'D', 10, 'errfncs', err, 'maxit', maxit);

%% AdaGrad
[~,infoAdaGrad_r] = AdaGradSDP( C, b, Aop, Atop, alpha, r, 'D', 10, 'errfncs', err, 'maxit', maxit);

[~,infoAdaGrad_3r] = AdaGradSDP( C, b, Aop, Atop, alpha, 3*r, 'D', 10, 'errfncs', err, 'maxit', maxit);

%% AdaNGD
[~,infoAdaNGD_r] = AdaNGDSDP( C, b, Aop, Atop, alpha, r, 1.1, 'D', 100, 'errfncs', err, 'maxit', maxit);

[~,infoAdaNGD_3r] = AdaNGDSDP( C, b, Aop, Atop, alpha, 3*r, 1.1, 'D', 100, 'errfncs', err, 'maxit', maxit);

%% Save Results
if ~exist('./results','dir'), mkdir results; end
save('./results/Test2_MatrixCompletion_Results.mat','info*');

%% Last edit: Alp Yurtsever - November 6, 2019
