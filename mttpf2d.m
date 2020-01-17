function [Xr Pr]=mttpf2d(Xr, Pr,  Zk, dt, calib)

global k

m=7;
n=10;
N=500;
ck=1;
cid=2;
ci=2;
gt=calib.center(1)/2;
da_type=0;

%%%%
if isempty(Xr)
    Xr=zeros(1,n+ci);
    Pr=zeros(1,n*N+ci);
end

%%% for adding a point
dist=[1 1 5 5 10 0 0 0 0 0]'; % units in [cm cm cm/s cm/s pixels]
if size(Xr,1)==1 && Xr(1)
    eta1=randn(n,N);
    P0=Xr(ci+1:n+ci)'*ones(1,N) + (dist*ones(1,N)).*eta1;
    Pr=[k Xr(2), P0(:)'];
    wts=normpdf(eta1(1,:), 0, 1);
end


if size(Zk,1) ~=m
    Zk=[Zk; zeros(m-size(Zk,1),size(Zk,2))];
end

if size(Xr,2) > n+2
    Xr=Xr(:,1:n+2);
end

% get the number of records in the current time step
curr=(Xr(:,ck)==k-1 & Xr(:,ck)~=0);
Pk=Pr(curr,ci+1:N*n+2);
P=[]; wts=[];
for ii=1:sum(curr)
    P(n*(ii-1)+1:n*ii,1:N)=reshape(Pk(ii,:),n,N);
    wts(ii,1:N)=ones(1,N);
end

nt=sum(curr);
ids=unique(Xr(:,cid));
if isempty(ids)
    idnew=1;
else
    idnew=max(ids)+1;
end


X0=[0, 0, 0, 0, 20, 1 zeros(1,4)];
P0=X0'*ones(1,N);

nZk=size(Zk,2);

% --- filter f_update (measurement model)
H=zeros(m,n); 

H(1,1)=1/calib.pix2cm(1); H(2,2)=1/calib.pix2cm(2); H(1:2,6)=calib.center;
H([3 4 5 6 7],[5 7 8 9 10])=eye(5); Hk=H;
eta=diag(1./calib.pix2cm); % note how we put the error in terms of the number of pixels per cm
Hinv=@(z) [(z(1:2)-calib.center).*(ones(2,1).*calib.pix2cm'); 0; 0; z(3); 1; z(4:7)]; % converts from Z to X
f_update=@(k, as1, P,  wts, Zk) pfUpdateMT(k, as1, P, wts, Zk, Hk, eta, n);


% --- f_predict (motion model)
Fk=eye(n); Fk(1,3)=dt; Fk(2,4)=dt; 
f_predict=@(k, P) pfPredictMT(k, P, n, Fk, dist);

initialize=@(k, Zk, idnew) init(k, Zk, Hinv, P0, n, ci, dist, idnew);

% --- associate
switch da_type
    case {0, 'gnn'}
        da1=@(D) gnnda(D, gt); 
    case {1, 'gnn+nnda'}
        da1=@(D) gnnda_nnda(D, gt); 
    case {3, 'nnda'}
        da1=@(D) nnda(D);
end


zhat=@(k, r, Xh) Hk*pfout(P(r, :),[], 1);

gvfun=@(Zh, Zk) norm(Zh-Zk);
costmatrix=@(Zh, Zk) costmatrix2d(Zh, Zk, gvfun);

Zh=zeros(m,nt);

% f_predict
if k>1
    P=f_predict(k, P);
end
for ii=1:nt
    [r c]=getind(n, 1, ii, 1:n, n);
    if P(r(1),1)
        Zh(:,ii)=zhat(1, r, P);
    end
end
Zh=Zh(:,Zh(1,:)~=0);

if ~isempty(Zk)
    D = costmatrix(Zh, Zk);

    % associate
    as1=da1(D);

    % f_update
    [P wts]=f_update(k, as1, P, wts, Zk);
else % if there is no measurement then terminate the target
    as1=0;
    Zk=zeros(m,1);
    % f_update
    [P wts]=f_update(k, as1, P, wts, Zk);    
end
% if sum(as1==0), keyboard; end
% initialize
if ~isempty(Zk)
Zk=Zk(:,~ismember(1:nZk,as1));
end
[Xi Pi]=initialize(k, Zk, idnew);


%%%%%%
% Xh=pfout(P, wts, 2);
Xh=pfout(P, wts, 1);
Xru(:,ci+1:n+ci)=reshape(Xh, n, sum(curr))';
Xru(:,ck)=k;
Xru(:,cid)=Xr(curr,cid);


for ii=1:sum(curr)
    P1=P(n*(ii-1)+1:n*ii,1:N);
    Prt(ii,:)=P1(:);
end
Pru=[];
if sum(curr)
    Pru(:,ci+1:n*N+ci)=Prt;
    Pru(:,ck)=k;
    Pru(:,cid)=Pr(curr,cid);
end


Xr=[Xr; Xru; Xi];
Pr=[Pr; Pru; Pi];

%%% delete all records that are five steps ago, current step, future, and
%%% not updated
delk=Xr(:,ck)<=k-5 | Xr(:,ck)>k | Xr(:,3)==0;
Xr(delk,:)=[];
Pr(delk,:)=[];


function D = costmatrix2d(Zh, Zk, gvfun)

nZ=size(Zk,2);
nT=size(Zh,2);
D=ones(nZ,nT)*1000000;
for tt=1:nT
    for zz=1:nZ
        D(zz,tt)=gvfun(Zh(1:2,tt), Zk(1:2,zz)); % position only
%         D(zz,tt)=gvfun(Zh(:,tt), Zk(:,zz), S(:,:,tt)); % position and size
    end
end

function [p wts] = pfUpdate(p, wts, Z, H, eta)
% function P = pfUpdate(P, Z, H, eta)

for jj=1:size(p,2)
    X=H*p(:,jj);

    % assuming independence along x and y dimensions
    wts(jj)=normpdf(X(1,:), Z(1), eta(1,1)).*...
         normpdf(X(2,:), Z(2), eta(2,2));
end

% for the rest just copy
p([5 7 8 9 10],:)=Z(3:7)*ones(1,size(p,2));
 
wts=wts/sum(wts);
neff=1/sum(wts.^2);
if (neff <=  size(p,2)/2)
    p=p(:,resample(wts));
end 

function [P wts]=pfUpdateMT(k, as1, P, wts, Zk, Hk, eta, n)

k=1;

% update
for tt=1:size(P,1)/n
    [r c]=getind(n, k, tt, 1:n, n);
    if numel(as1) >= tt
        if as1(tt)
            [P(r,:), wts(tt,:)]=pfUpdate(P(r,:), wts(tt,:), Zk(:,as1(tt)), Hk, eta);
        else
            P(r,:)=0;
            wts(tt,:)=1;
        end
    end
end

function p = pfPredict(p, F, dist, n)
%function [Xh, P] = kalmanPredict(Xh, P, Klmn)
%
% F is nxn transition matrix
% Q is nxn disturbance matrix

for jj=1:size(p,2)
    p(:,jj)=F*p(:,jj) + dist.*randn(n,1);
end

function P= pfPredictMT(k, P, n, F, dist)

k=1;
for ii=1:size(P,1)/n
    [r c]=getind(n, k, ii, 1:n, n);
    if P(r(1),1)
        P(r,:)=pfPredict(P(r,:), F, dist, n);
    end
end

function [Xi Pi]=init(k, Zk, Hinv, P0, n, ci, dist, idnew)
nZk=size(Zk,2);
Pi=zeros(nZk,numel(P0)+ci);
Xi=zeros(nZk,n+ci);

% initialize
if ~isempty(Zk)
    for zz=1:nZk
        P1=P0;
        for jj=1:size(P1,2)
            P1(:,jj)=Hinv(Zk(:,zz)) + dist.*randn(n,1);
        end
        Xi(zz,:)=[k, idnew, pfout(P1,[],1)'];
        Xi(zz,5:6)=0; % we want initial velocity to be zero
        Pi(zz,:)=[k, idnew, P1(:)'];
        idnew=idnew+1;
    end
end


function as1 = gnnda(D, thresh)

as=munkres(D);
idx=find(as>0);

% remove all associations beyond a certain threshold
as(idx(D(idx)>thresh))=0; % as is number of targets with each having a measurement number

[val, as1]=max(as,[],1);
as1(val==0)=0;

function as1 = nnda(D)

[jnk as1]=min(D,[],1);


function as1 = gnnda_nnda(D, thresh)


as=munkres(D);
idx=find(as>0);

% remove all associations beyond a certain threshold
as(idx(D(idx)>thresh))=0; % as is number of targets with each having a measurement number

[val, as1]=max(as,[],1);
as1(val==0)=0;

%%%% nearest neighbor
for ii=find(as1==0)
    [vv zz]=min(D(:,ii));
    if vv<thresh/2
        as1(ii)=zz;
    end
end



function estx = pfout(x, wts, flag)

%{
This function will return an estimate of a distribution, based on what flag is passed.
So, for e.g. if flag passed is 1 and x is a matrix, the function will compute the mean along
the longer dimension
if flag is 2 then the function will compute the max along the longer dimension... and so on.

ok. so there's a check now that the longer dimension has to be columns!!
%}

% if all wts are same
% if ~max(wts), flag=1; end
if ~isempty(x)
    switch flag
        case 1
            estx = mean(x,2);
        case 2
            n=size(x,1)/size(wts,1);
            for jj=1:size(wts,1)
                [val idx] = max(wts(jj,:));
                estx(n*(jj-1)+1:n*jj,1) = x(n*(jj-1)+1:n*jj,idx);
            end
        case 3
            estx = mode(x,2);
        otherwise
            error('brrrr!!');
    end
else
    estx=[];
end

function S=zhcov(r,P, Hk)
zh=zeros(size(Hk,1),size(P,2));
for jj=1:size(P,2)
    zh(:,jj)=Hk*P(r,jj);
end
S=cov(zh');