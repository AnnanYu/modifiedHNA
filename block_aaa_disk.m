function [R,rmse,opts,A,B,C,t,H,G,maxerr] = block_aaa_disk(F,pts,opts)
%BLOCK_AAA   Block-AAA algorithm for discrete rational LS approximation.
%
% Block generalization of the AAA algorithm presented in 
% [Nakatsukasa/Sete/Trefethen, SIAM J. Sci. Comput. 40 (3), 2018]
% to m-by-n functions. Block-AAA is based on a generalized barycentric
% formula with matrix-valued weights, that is,
%
%    R(z) = (sum_i W_i/(z-z_i))\(sum_i W_i F(z_i)/(z-z_i))
%
% where the W_i are m-by-m matrices and z_i are scalar support points. 
% Note that the F(z_i) are m-by-n matrices and hence the numerator 
% effectively uses m-by-n matrices for interpolation (if W_i is nonsing.)
%
% Block-AAA is called as [R,rmse,out] = block_aaa(F,pts,opts)
%
% Inputs:  F    -- function handle @(z) to the m-by-n function, 
%                  or a cell array of length(pts) m-by-n matrices 
%          pts  -- vector of distinct sampling points in the complex plane
%          opts -- structure of additional options (default values):
%                  opts.tol    : target root mean squared error (1e-12)
%                  opts.maxit  : maximal number of iterations 
%                  opts.return : which baryfun to return, (best), last, all.
%
% Returns: R    -- rational interpolant as a BARYFUN object
%          rmse -- root mean squared error for each block-AAA iteration
%          out  -- additional outputs.
%
% Website: see https://github.com/nla-group/block_aaa for more details  
%

% handle input options
if nargin < 3, opts = struct(); end
if ~isfield(opts,'tol'), opts.tol = 1e-12; end
if ~isfield(opts,'reg'), opts.reg = 1e-12; end
if ~isfield(opts,'maxit'), opts.maxit = 50; end
if ~isfield(opts,'svd'), opts.svd = 1; end % use SVD, otherwise EIG (not recomm.)
if ~isfield(opts,'chol'), opts.chol = 0; end % use QR+Cholesky update, otherwise full SVD
if ~isfield(opts,'return'), opts.return = 'best'; end % can be 'best', 'last', 'all'

pts = pts(:);
npts = length(pts);

% if F is function handle, evaluate it at all the pts
if iscell(F)
    FF = F;
else
    FF = cell(npts,1);
    for i = 1:npts
        FF{i} = F(pts(i));
    end
end
[m,n] = size(FF{1});

t = 0;

% compute norm ||F(z)|| at all pts
err = 0*pts;
for i = 1:npts
    err(i) = norm(FF{i},'fro').^2;
end
% find max and assign first support point there
[~,ind] = max(err);

lamind = 1:npts;              % indices of approximation points
lamind(lamind==ind) = [];
zk_ind = ind;
lam = pts(lamind);
%lam = pts;
zk = pts(zk_ind);

% M is the block Loewner matrix
M = zeros(0,(npts-1)*n);
%M = zeros(0,npts*n);

if strcmp(opts.return,'all'), R = cell(1); else R = 0; end

for it = 1:opts.maxit+1
    
    ell = length(zk)-1; % degree
    p = npts-ell-1;
    %p = npts;
    B = zeros(m,p*n);
    for i = 1:p % add new row to M
            M(ell*m+1:(ell+1)*m,(i-1)*n+1:i*n) =...
                (FF{lamind(i)} - FF{zk_ind(ell+1)})/(lam(i) - zk(ell+1));
        B(:,(i-1)*n+1:i*n) = -FF{(i)};
    end
%     [U,S,V] = svd(M'); tol = 1e-5;
%     ind = find(diag(S) > tol);
%     U = U(:,ind); S = S(ind,ind); V = V(:,ind);
%     EF = (V * diag(diag(S).^-1) * U' * B')';
    tic
    EF = (M'\B')';
    t = t + toc;
%     [U,S,V] = svd(M'); tol = 1e-6;
%     ind = find(diag(S) > tol);
%     U = U(:,ind); S = S(ind,ind); V = V(:,ind);
%     M = U*S*V'; M = M';
%    EF = lsqminnorm(M',B'); EF = EF';
    %EF = (pinv(M')*B')';
    %pen = opts.reg; [U,S,V] = svd(M');
    %EF = (V * diag((diag(S).^2 + (pen/it)).^(-1)) * S' * (U' * B'))';
    
    
    % construct function handle to barycentric evaluation
    Ck = cell(ell+1,1); Dk = cell(ell+1,1);
    for j = 1:ell+1
        Dk{j} = EF(:,(j-1)*m+1:j*m);
        Ck{j} = Dk{j}*FF{zk_ind(j)};
    end
    % we don't use baryfun class here for performance
    R = @(z) eval_bary(z,zk,Ck,Dk); 
    
    % evaluate RMSE
    err = 0*pts;
    for i = 1:npts
        if ~ismember(i,zk_ind)
            err(i) = norm(R(pts(i)) - FF{i},'fro');%.^2;
        end
    end
    
    indnan = find(not(isfinite(err)));
    err(indnan) = 0;
    
    rmse(it) = sqrt(sum(err)/npts);
    maxerr(it) = max(err);
    
    if rmse(it) < opts.tol || it == opts.maxit
        break
    end
    
    % not done: find max error and assign next support point
    [~,ind] = max(err);
    ind2 = find(lamind==ind(1));
    lamind(ind2) = [];
    zk_ind = [zk_ind , ind];
    lam = pts(lamind);
    zk = pts(zk_ind);

    M(:,1+(ind2-1)*n:ind2*n) = [];
    
end % for it
opts.zk = zk;

%t = toc;

tic

B  = []; C = [];
for i = 1:it
    B = [B;FF{zk_ind(i)}]; C = [C,Dk{i}];
end
A = kron(diag(zk),eye(m)) - kron(ones(it,1),eye(m))*C;
%[A,B,C,D] = plane2disk(A,B,C,zeros(m,n));
[X,Lam] = eig(A);
F1 = @(z) C * (((z)*eye(length(A)) - A) \ B);
F2 = @(z) R(z);
ind_stable = find(abs(diag(Lam)) < 1-1e-8);
ind_unstable = find(abs(diag(Lam)) >= 1); XinvB = X \ B;
%ind_stable = find(abs(diag(Lam)) < 1-1e-8);
%ind_unstable = find(abs(diag(Lam)) >= 1); XinvB = X \ B;
H = @(z) C * X(:,ind_unstable) * ...
    diag(diag((z*eye(length(ind_unstable)) - Lam(ind_unstable,ind_unstable))).^-1) * ...
    XinvB(ind_unstable,:);
A = Lam(ind_stable,ind_stable); C = C * X(:,ind_stable); B = XinvB(ind_stable,:);

G = @(z) C * ((z*eye(length(A)) - A) \ B);
F3 = @(z)G(z) + H(z);

t = t + toc;


end % block_aaa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = eval_bary(z,zk,Ck,Dk)
%EVAL_BARY   Evaluate matrix-valued barycentric formula
%   
%      R(z) = inv(sum_k Dk/(z-zk))*(sum_k Ck/(z-zk))
%
% where z is the scalar evaluation point, zk is a vector of the distinct
% support points, and Ck and Dk are cell arrays with the corresponding
% numerator and denominator coefficients, respectively. 

N = zeros(size(Ck{1})); % numerator
D = zeros(size(Dk{1})); % denominator

[val,ind] = min(abs(z-zk)); % evaluation at support point
if val < 1e1*eps
    R = (Dk{ind})\Ck{ind};
    return
end

for j = 1:length(zk)
    N = N + Ck{j}/(z-zk(j));
    D = D + Dk{j}/(z-zk(j));
end
R = (eye(length(D)) + D)\N;

end
