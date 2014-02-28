function [Xr Pr]=mtt2d(Xr, Pr,  Zk, dt, calib, da_type)

global k

m=3;
n=6;
ck=1;
cid=2;
ci=2;
gt=9;

%%%%
curr=(Xr(:,ck)==k-1 & Xr(:,ck)~=0);

Xk=Xr(curr,ci+1:n+2)'; % this is just the current time-step
Xh=Xk(:);
Pk=Pr(curr,ci+1:n^2+2);
P=[];
for ii=1:sum(curr)
    P(n*(ii-1)+1:n*ii,1:n)=reshape(Pk(ii,:),n,n);
end

nt=sum(curr);
ids=unique(Xr(:,cid));
if isempty(ids)
    idnew=1;
else
    idnew=max(ids)+1;
end


X0=[1, 1, 0, 0, 1, 1];
P0=eye(n);

nZk=size(Zk,2);

% --- filter f_update (measurement model)
H=zeros(m,n); 

H(1,1)=1/calib.pix2cm(1); H(2,2)=1/calib.pix2cm(2); H(1:2,n)=calib.center;
H(3,5)=1; Hk=H;
R=eye(m); R(1:2,1:2)=diag(1./calib.pix2cm)*2; % note how we put the error in terms of the number of pixels per cm
R(m,m)=20;% measurement noise sqrt of these values is measurement noise
Hinv=@(z) [(z(1:2)-calib.center).*(ones(2,1).*calib.pix2cm'); 0; 0; z(3); 1]; % converts from Z to X

f_update=@(k, as1, Xh, Zk) UpdateMT(k, as1, Xh, Zk, Hinv);
initialize=@(k, Zk, Xh, P, idnew) init(k, Zk, Hinv, X0, n, ci, idnew);

% --- f_predict (motion model)
[F, Fk]=motion_model(dt, n);
f_predict=@(k, Xh) PredictMT(k, Xh, F, n);

% --- associate
switch da_type
    case {0, 'gnn'}
        da1=@(D) gnnda(D, gt); 
    case {1, 'gnn+nnda'}
        da1=@(D) gnnda_nnda(D, gt); 
end

% --- measurement model
zhat=@(k, r, Xh) Hk*Xh(r, k);

% --- measurement covariance
gvfun=@(Zh, Zk) norm(Zh-Zk);
costmatrix=@(Zh, Zk) costmatrix2d(Zh, Zk, gvfun);

Zh=zeros(m,nt);
S=zeros(m,m,nt);

% prediction / motion model
if k>1
    [Xh, P]=f_predict(k, Xh);
end
for ii=1:nt
    [r c]=getind(n, 1, ii, 1:n, n);
    if Xh(r(1),1)
        Zh(:,ii)=zhat(1, r, Xh);
    end
end
Zh=Zh(:,Zh(1,:)~=0);

% matching measurements to targets
if ~isempty(Zk)
    D = costmatrix(Zh, Zk);

    % associate 
    as1=da1(D);

    % f_update
    [Xh P]=f_update(k, as1, Xh, Zk);
else % if there is no measurement then terminate the target
    as1=0;
    Zk=zeros(m,1);
    % f_update
    [Xh P]=f_update(k, as1, Xh, Zk);
end

% initialize
if ~isempty(Zk)
	Zk=Zk(:,~ismember(1:nZk,as1));
end
[Xi]=initialize(k, Zk, Xh,idnew);


% put it back in the output form of X
Xru(:,ci+1:n+ci)=reshape(Xh, n, sum(curr))';
Xru(:,ck)=k;
Xru(:,cid)=Xr(curr,cid);


% delete all records that are five steps ago, current step, future, and
% not updated
delk=Xr(:,ck)<=k-5 | Xr(:,ck)>k | Xr(:,3)==0;
Xr(delk,:)=[];


function D = costmatrix2d(Zh, Zk, S, gvfun)

nZ=size(Zk,2);
nT=size(Zh,2);
D=ones(nZ,nT)*1000000;
for tt=1:nT
    for zz=1:nZ
        D(zz,tt)=gvfun(Zh(1:2,tt), Zk(1:2,zz), S(1:2,1:2,tt)); % position only
%         D(zz,tt)=gvfun(Zh(:,tt), Zk(:,zz), S(:,:,tt)); % position and size
    end
end

function [F, Fk]= motion_model(dt,n)

Fk=eye(n); Fk(1,3)=dt; Fk(2,4)=dt; 
F=Fk; 

function Xh = Update(Xh, Z, Hinv)

Xh=Hinv(Z);

function Xh=UpdateMT(k, as1, Xh, Zk, Hinv)

k=1;

% update
for tt=1:size(Xh,1)/n
    [r c]=getind(n, k, tt, 1:n, n);
    if numel(as1) >= tt
        if as1(tt)
            Xh(r,k)=Update(Xh(r,k), Zk(:,as1(tt)),  Hinv);
        else
            Xh(r,k)=0;
            P(r,c)=0;
        end
    end
end

function [Xh] = Predict(Xh,F)
%function [Xh, P] = Predict(Xh, P, Klmn)
%
% F is nxn transition matrix
% Q is nxn disturbance matrix


% state
Xh= F*Xh;


function [Xh]= PredictMT(k, Xh,F, n)

k=1;
for ii=1:size(Xh,1)/n
    [r c]=getind(n, k, ii, 1:n, n);
    if Xh(r(1),1)
        [Xh(r,k)]=Predict(Xh(r,1), F, Q);
    end
end

function [Xi Pi]=init(k, Zk, Hinv, X0, n, ci, idnew)
nZk=size(Zk,2);
Xi=zeros(nZk,n+ci);
Pi=zeros(nZk,n*n+ci);

% initialize
if ~isempty(Zk)
    for zz=1:nZk
        Xi(zz,:)=[k, idnew,  X0]; 
        Xi(zz, ci+1:n+ci)=Hinv(Zk(:,zz));
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


