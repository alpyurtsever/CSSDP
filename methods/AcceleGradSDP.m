function [ y_bar, out ] = AcceleGradSDP( C, b, Aop, Atop, alpha, r, varargin )
%ACCELEGRADSDP Implements AcceleGrad [LYC2018], for solving the model dual 
%SDP problem template in [DYTCU2019]. 
%
%ACCELEGRADSDP requires at least 5 inputs: 
% - (sparse) cost matrix C
% - measurement (constraint) vector b
% - function handle Aop for constraint matrix A: 
%   Aop(U,S) should compute A*(U*S*U') efficiently, without constructing
%   (U*S*U') in the ambient dimensions.
% - function handle Atop for constraint matrix A: 
%   Atop(y,z) should compute (A'*y)*z efficiently, without constructing
%   (A'*y) in the ambient dimensions.
% - scalar r for the reconstruction rank
%
%ACCELEGRADSDP also has optional property inputs. These properties should  
%be input as a tuple {'property-tag',value} after 5th input. For example 
%AcceleGradSDP(C, b, Aop, Atop, alpha, r, 'maxit', 100, 'eigstol', 1e-3)
%The list of properties and the default values are as follows:
%
% - 'x0' (default: zeros)
%   Sets the initial estimate.
%
% - 'maxit' (default: 1000)
%   Sets the number of iterations to run.
%
% - 'D' (default: 1)
%   Parameter for the diameter, effects the step-size.
%
% - 'doubling' (default: false)
%   Enables the heuristic doubling rule for 'D' if set to true.
%
% - 'G' (default: 0)
%   Parameter for the Lipschitz constant, effects the step-size. 
%    
% - 'eigsTol' (defalut: 1e-6)
%   Sets the tolerance parameter for the eigenvector subsolver for 
%   subgradient oracle. 
%
% - 'errfncs' (default: empty)
%   Cell for function handles to compute errors. Odd inputs entries of the
%   cell sets a name for the error function, and the even entries are the
%   function handles to compute errors. For example, 
%   err{1} = 'NormX';
%   err{2} = @(y,S,V) norm(V*S*V','fro');
%   err{3} = 'Normy';
%   err{4} = @(y,S,V) norm(y,2);
%   If we input {'errfncs',err} this will compute the Frobenius norm of the 
%   primal estimate X in the factored form V*S*V', and the 2-norm of the 
%   dual estimate y. These errors are written to out.y_bar.NormX and 
%   out.y_bar.Normy. The function handle should always have 3 inputs. The 
%   first input is to compute errors from the dual estimate, the 2nd and  
%   3rd inputs are to compute errors from the primal estimate. 
%
% - 'saveit' (default: a combination of linear and logarithmically spaced 
%   iterations - see lines 75-78)
%   We run the Complementary Slackness Primal Recovery only at some
%   iterations to avoid the cost. 'saveit' determines at which iterations
%   we will construct a primal estimate. Errors and convergence measures
%   are computed only at these iterations.
%
% - 'saveName' (default: empty)
%   Path to save intermediate results. For example, if you set 'saveName'
%   to 'results/inter.mat' the results will be saved to this location/file 
%   at iterations set by 'saveit'
%
%ACCELEGRADSDP has two outputs
%
% - y_bar : the dual estimate
%   We output the averaged sequence with convergence guarantees.
%   Heuristically, y_t also converges in most scenarios, sometimes
%   significantly faster than y_bar. (You need to modify the code to output
%   y_t or to compute errors from y_t).
%
% - out : Struct that keeps the algorithmic and convergence information.
%   "out.info" - keeps the DIMENSION, reconstruction RANK (input r), and a
%   flag SUBSOLVERCONVERGES. SUBSOLVERCONVERGES is true if eigs converged
%   at all iterations, and false if it fails to achieve the target accuracy
%   at any iteration. 
%
%   "out.params" - keeps the information about the parameter setting : 
%   ALPHA, ACCELEGRAD_D, ACCELEGRAD_G, and FLAG_DOUBLING.
%   See (input 'alpha' and properties 'D','G','doubling').
%
%   "out.iteration" - keeps a vector indicating the iterations at which we 
%   construct the primal estimate and measure the error. 
%   See (property 'saveit').
%
%   "out.time" - keeps timing information. 
%   <dualSolver> shows the amount of time used to solve the dual.
%   <compslackRec> shows the amount of time used to construct the primal.
%
%   "out.y_bar" - keeps error information computed from the averaged dual 
%   sequence and the primal recovery from it. 
%   <dualObj> (constrained) dual objective "y'*b"
%   <dualFeas> (constrained) dual infeasibility from the dual slack matrix
%   <dualTot> complete (unconstrained) dual objective that we actually use
%   <primalObj> primal objective value "trace(C'*X)"
%   <primalFeas> primal infeasibility "norm(AX-b)"
%   <optional> you can add any other error using the 'errfncs' property.
%   
% [LYC2018] K. Levy, A. Yurtsever, V. Cevher,
% "Online Adaptive Methods, Universality and Acceleration" 
% Advances in Neural Information Processing Systems 31, (NeurIPS 2018).
%
% [DYTCU2019] L. Ding, A. Yurtsever, V. Cevher, J.A. Tropp, M. Udell,
% "An Optimal-Storage Approach to Semidefinite Programming using Approximate Complementarity"
% arXiv:1902.03373, 2019.

optsEigs.isreal = isreal(C);      % We are working on real number field.
optsEigs.issym  = issymmetric(C); % Gradient is symmetric.
optsEigs.tol    = 1e-6;           % Stopping tolerance for eigs.
optsEigs.p      = 5;
optsEigs.maxit  = 500;

sName = [];
err = {};
x0 = zeros(size(b));
D = 1;
G = 0;
FLAG_DOUBLING = 0;
maxit = 1e3;
if ~isempty(varargin)
    for tt = 1:2:length(varargin)
        switch lower(varargin{tt})
            case 'errfncs'
                err = varargin{tt+1};
            case 'savename'
                sName = varargin{tt+1};
            case 'x0'
                x0 = varargin{tt+1};
            case 'doubling'
                FLAG_DOUBLING = varargin{tt+1};
            case 'd'
                D = varargin{tt+1};
            case 'g'
                G = varargin{tt+1};
            case 'maxit'
                maxit = varargin{tt+1};
            case 'saveit'
                saveHist = varargin{tt+1};
            case 'eigstol'
                optsEigs.tol = varargin{tt+1};
        end
    end
end
saveHistLin = unique(round(linspace(1,maxit,30)));
saveHistLog = unique(round(1.25.^(0:1:100)));
saveHistLog(saveHistLog>maxit) = [];
saveHist = unique(sort([0,1,saveHistLin,saveHistLog,maxit]));

n = size(C,1);
x_t = x0;
y_t = x0;
z_t = x0;
y_bar = x0;

S1 = 0;
S2 = 0;

[out, errNames, ptr] = createErrStructs();
printHeader(errNames);

out.info.SUBSOLVERCONVERGES = true;
out.info.RANK = r;
out.info.DIMENSION = n;
out.params.ACCELEGRAD_D = D;
out.params.ACCELEGRAD_G = G;
out.params.ALPHA = alpha;
out.params.FLAG_DOUBLING = FLAG_DOUBLING;

% Timer for dual solver
TimeDual = 0;

% Loop
for t = 0:maxit
    
    % Start timer
    clkDual = tic;
    
    alpha_t = max(0.25*(t+1),1);
    % Set
    tau_t = 1/alpha_t;
    
    % Update
    x_t = tau_t*z_t + (1-tau_t)*y_t;
    
    % Compute the subgradient
    eigsArg = @(x) C*x - Atop(x_t,x);
    for tt = 1:20
        [u_t,sigma_t,eigsFlag] = eigs(eigsArg, n, 1, 'sa', optsEigs);
        if ~eigsFlag, break; else
            optsEigs.p = optsEigs.p + 1;
            if optsEigs.p > n
                warning('Subsolver is not converging!');
                out.info.SUBSOLVERCONVERGES = false;
                optsEigs.p = n;
                break;
            end
        end
    end
    if tt == 20
        warning('Subsolver is not converging!');
        out.info.SUBSOLVERCONVERGES = false;
    end
    
    if sigma_t < 0
        g = -b + alpha*(Aop(u_t,1));
    else
        g = -b;
    end
    
    ng = norm(g);
    
    S1 = S1 + (alpha_t*ng)^2;
    
    eta_t = D / sqrt(2*(G^2+S1));
    
    z_t = z_t - alpha_t*eta_t*g;
    
    y_t = x_t - eta_t*g;
    if FLAG_DOUBLING
        ny = norm(y_t);
        if ny > D
            fprintf('| %7d |', t );
            fprintf(' *** increase diameter from %d ',D);
            D = 2*D;
            fprintf('to %d \n',D);
        end
    end
    
    % Averaged
    S2 = S2 + alpha_t;
    weight = alpha_t/S2;
    y_bar = (1 - weight)*y_bar + weight*y_t;
    
    % Stop timer
    TimeDual = TimeDual + toc(clkDual);
    
    if any(t==saveHist)
        
        % recover primal
        clkRec = tic;
        
        eigsArg = @(x) C*x - Atop(y_bar,x);
        for ttt = 1:20
            [u_t,sigma_t,eigsFlag] = eigs(eigsArg, n, 1, 'sa', optsEigs);
            if ~eigsFlag, break; else
                optsEigs.p = optsEigs.p + 1;
                if optsEigs.p > n
                    warning('Subsolver is not converging!');
                    out.info.SUBSOLVERCONVERGES = false;
                    optsEigs.p = n;
                    break;
                end
            end
        end
        if ttt == 20
            warning('Subsolver is not converging!');
            out.info.SUBSOLVERCONVERGES = false;
        end
        
        [Vy,Sy,Dy] = recoverPrimal(@(x) C*x - Atop(y_bar,x), n, r,...
            optsEigs, Aop, b); %#ok
        
        TimeCompslackRec = toc(clkRec);
        
        updateErrStructs();
        printError();
        saveFile();
        
    end
end

%% End of main algorithm
% Nested functions


    function [out, errNames, ptr] = createErrStructs()
        out.time.dualSolver = nan(numel(saveHist),1);
        out.time.compslackRec = nan(numel(saveHist),1);
        
        out.y_bar.dualObj = nan(numel(saveHist),1);
        out.y_bar.dualFeas = nan(numel(saveHist),1);
        out.y_bar.dualTot = nan(numel(saveHist),1);
        out.y_bar.primalObj = nan(numel(saveHist),length(r));
        out.y_bar.primalFeas = nan(numel(saveHist),length(r));
        
        for eIt = 1:2:length(err)
            out.y_bar.(err{eIt}) = nan(numel(saveHist),length(r));
        end
        
        out.iteration = nan(numel(saveHist),1);
        errNames = fieldnames(out.y_bar);
        ptr = 0;
    end

    function updateErrStructs()
        ptr = ptr+1;
        out.iteration(ptr,1) = t;
        
        out.time.dualSolver(ptr,1) = TimeDual;
        out.time.compslackRec(ptr,1) = TimeCompslackRec;
        
        out.y_bar.dualObj(ptr,1) = y_bar'*b;
        out.y_bar.dualFeas(ptr,1) = abs(min(sigma_t,0));
        out.y_bar.dualTot(ptr,1) = y_bar'*b + alpha*min(sigma_t,0);
        out.y_bar.primalObj(ptr,1) = sum(sum((C*Vy).*(Vy*Sy)));
        out.y_bar.primalFeas(ptr,1) = norm(Aop(Vy,Sy) - b,'fro');
        
        for eIt = 1:2:length(err)
            out.y_bar.(err{eIt})(ptr,1) = err{eIt+1}(y_bar,Sy,Vy);
        end
    end

    function printError()
        fprintf('| %7d |', t );
        for pIt = 1:length(errNames)
            if ~isnan(out.y_bar.(errNames{pIt})(ptr,1))
                fprintf(' % 5.3e |', out.y_bar.(errNames{pIt})(ptr,1) );
            else
                fprintf(' % 10s |', 'NaN' );
            end
        end
        fprintf('\n');
    end

    function saveFile()
        if ~isempty(sName)
            save(sName,'out','-v7.3');
        end
    end

end

function [V,S,D] = recoverPrimal(eigsOp, n, r, optsEigs, Aop, b) %#ok
optsEigs.p = min(2*r,n);
for tt = 1:20
    [V,D,eigsFlag] = eigs(eigsOp ,n ,r ,'sa', optsEigs);
    if ~eigsFlag, break; else
        optsEigs.p = optsEigs.p + r;
        if optsEigs.p > n
            warning('Primal recovery not accurate!');
            optsEigs.p = n;
            break;
        end
    end
end
cvx_begin
cvx_quiet true
variable S(r,r) semidefinite
minimize  norm(Aop(V,S)-b,'fro')
cvx_end
end

function printHeader(errNamesPrint)
for pIt = 1:length(errNamesPrint)
    if length(errNamesPrint{pIt}) > 10
        errNamesPrint{pIt} = errNamesPrint{pIt}(1:10);
    end
end
fprintf('|    iter | ');
for pIt = 1:length(errNamesPrint)
    fprintf('%10s', errNamesPrint{pIt});
    fprintf(' | ');
end
fprintf('\n');
end

%% Last edit: Alp Yurtsever - November 6, 2019