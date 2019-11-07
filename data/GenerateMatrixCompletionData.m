%% Description
% This script draws a random data for the Matrix Completion experiments. 
% It saves the data as "MatrixCompletionData.mat". 
% We include the output of this script in the repository since it requires
% large storage resources to run it. Therefore, you do not need to run this
% script to reproduce our results in [DYTCU2019]. 
%
% [DYTCU2019] L. Ding, A. Yurtsever, V. Cevher, J.A. Tropp, M. Udell,
% "An Optimal-Storage Approach to Semidefinite Programming using Approximate Complementarity"
% arXiv:1902.03373, 2019.

%%% fix the random seed
clearvars;
rng(0,'twister');

%%% setup parameter
n_o = 1000;
m_o = 1500;
prob_o = 0.025;
n = m_o + n_o;                  % lifted dimensions
r = 5;                          % rank of X0

%%% Generating data
w = sign(randn(n_o,r));
v = sign(randn(m_o,r));
X = w*(v');
[u1,d1] = qr(w,0);
[u2,d2] = qr(v,0);
[u3,d3,v3] = svd(d1*d2');
U = u1*u3;
V = u2*v3;
D = d3;
nucNorm_X0 = sum(diag(d3));
O = rand(n_o,m_o);          % O: mask for observation probability
O = O <= prob_o;
Y = X.*O;                   % Y: Observed entries
[i,j,b]=find(Y);            % indices and values of observed entries
k = find(Y);                % indices and values of observed entries
m = length(b);

cvx_begin
cvx_solver sdpt3
cvx_precision best
variable Xcvx(n_o,m_o)
minimize norm_nuc(Xcvx)
subject to
Xcvx(k) == b;
cvx_end

% Lift the dimensions
[UU,DD,VV] = svd(Xcvx,'econ');
GTruth = [UU*DD*UU' Xcvx ; Xcvx' VV*DD*VV'];
Ylift = [zeros(n_o,n_o) Y; Y' zeros(m_o,m_o)];
[SubI,SubJ,Vals]=find(Ylift);
% Xlift = [U*D*U' X ; X' V*D*V'];
nucNorm_Xlift = 2*sum(diag(D));   % nuclear norm of X0
alpha = 1.1*nucNorm_Xlift;

save('./MatrixCompletionData.mat','SubI','SubJ','Vals','alpha','GTruth','r');

%% Last edit: Alp Yurtsever - November 6, 2019
