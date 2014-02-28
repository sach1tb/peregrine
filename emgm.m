function [label, model, llh] = emgm(X, init)
% Perform EM algorithm for fitting the Gaussian mixture model.
%   X: d x n data matrix
%   init: k (1 x 1) or label (1 x n, 1<=label(i)<=k) or center (d x k)
% Written by Michael Chen (sth4nth@gmail.com).
%% initialization
% fprintf('EM for Gaussian mixture: running ... \n');
R = initialization(X,init);
[jnk,label(1,:)] = max(R,[],2);
R = R(:,unique(label));

tol = 1e-6;
maxiter = 500;
llh = -inf(1,maxiter);
converged = false;
t = 1;
while ~converged && t < maxiter
    t = t+1;
    model = maximization(X,R);
    [R, llh(t)] = expectation(X,model);
    
    [jnk,label(1,:)] = max(R,[],2);
    idx = unique(label);   % non-empty components
    if size(R,2) ~= size(idx,2)
        R = R(:,idx);   % remove empty components
    else
        converged = llh(t)-llh(t-1) < tol*abs(llh(t));
    end

end
llh = llh(2:t);
if converged
%     fprintf('Converged in %d steps.\n',t-1);
else
%     fprintf('Not converged in %d steps.\n',maxiter);
end

function R = initialization(X, init)
[d,n] = size(X);
if isstruct(init)  % initialize with a model
    R  = expectation(X,init);
elseif length(init) == 1  % random initialization
    k = init;
    idx = randsample(n,k);
    m = X(:,idx);
    [jnk,label] = max(bsxfun(@minus,m'*X,sum(m.^2,1)'/2),[],1);
    while k ~= unique(label)
        idx = randsample(n,k);
        m = X(:,idx);
        [jnk,label] = max(bsxfun(@minus,m'*X,sum(m.^2,1)'/2),[],1);
    end
    R = full(sparse(1:n,label,1,n,k,n));
elseif size(init,1) == 1 && size(init,2) == n  % initialize with labels
    label = init;
    k = max(label);
    R = full(sparse(1:n,label,1,n,k,n));
elseif size(init,1) == d  %initialize with only centers
    k = size(init,2);
    m = init;
    [jnk,label] = max(bsxfun(@minus,m'*X,sum(m.^2,1)'/2),[],1);
    R = full(sparse(1:n,label,1,n,k,n));
else
    error('ERROR: init is not valid.');
end

function [R, llh] = expectation(X, model)
mu = model.mu;
Sigma = model.Sigma;
w = model.weight;

n = size(X,2);
k = size(mu,2);
logR = zeros(n,k);

for i = 1:k
    logR(:,i) = loggausspdf(X,mu(:,i),Sigma(:,:,i));
end
logR = bsxfun(@plus,logR,log(w));
T = logsumexp(logR,2);
llh = sum(T)/n; % loglikelihood
logR = bsxfun(@minus,logR,T);
R = exp(logR);


function model = maximization(X, R)
[d,n] = size(X);
k = size(R,2);

s = sum(R,1);
w = s/n;
mu = bsxfun(@times, X*R, 1./s);
Sigma = zeros(d,d,k);
for i = 1:k
    Xo = bsxfun(@minus,X,mu(:,i));
    Xo = bsxfun(@times,Xo,sqrt(R(:,i)'));
    Sigma(:,:,i) = Xo*Xo'/s(i);
    Sigma(:,:,i) = Sigma(:,:,i)+eye(d)*(1e-6); % add a prior for numerical stability
end

model.mu = mu;
model.Sigma = Sigma;
model.weight = w;

function y = loggausspdf(X, mu, Sigma)
d = size(X,1);
X = bsxfun(@minus,X,mu);
[R,p]= chol(Sigma);
if p ~= 0
    error('ERROR: Sigma is not PD.');
end
q = sum((R'\X).^2,1);  % quadratic term (M distance)
c = d*log(2*pi)+2*sum(log(diag(R)));   % normalization constant
y = -(c+q)/2;

function s = logsumexp(x, dim)
% Compute log(sum(exp(x),dim)) while avoiding numerical underflow.
%   By default dim = 1 (columns).
% Written by Michael Chen (sth4nth@gmail.com).
if nargin == 1, 
    % Determine which dimension sum will use
    dim = find(size(x)~=1,1);
    if isempty(dim), dim = 1; end
end

% subtract the largest in each column
y = max(x,[],dim);
x = bsxfun(@minus,x,y);
s = y + log(sum(exp(x),dim));
i = find(~isfinite(y));
if ~isempty(i)
    s(i) = y(i);
end