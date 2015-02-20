function [Xr Pr]=mttkf2ds(Xr, Pr,  Zk, dt, calib, da_type)
%function [Xr Pr]=mttkf2ds(Xr, Pr,  Zk, dt, calib, da_type)
%
% 2D kalman filter
% Xr is the [* x 6+2] row matrix for the current time step from the csv file
% Pr is the [* x 36+2] corresponding rows that will make the covariance matrix
% Zk is the [3 x *] matrix measurement for the current time step
% dt is the time between frames (e.g. 1/30 seconds)
% calib is the calibration matrix for the camera
% da_type is for data association, set it to 0

global k
global trklen

trklen=5;

m=3+4; 
n=6+4;
ck=1;
cid=2;
ci=2;
gt=16;

%%%%
if isempty(Xr)
    Xr=zeros(1,n+ci);
    Pr=zeros(1,n^2+ci);
end

%%% for adding a point
if size(Xr,1)==1 && Xr(1)
    P0=eye(n)*2; P0(3,3)=10; P0(4,4)=10; P0(5:n, 5:n)=eye(6);
    Pr=[k Xr(2), P0(:)'];
end

if size(Xr,2) > n+2
    Xr=Xr(:,1:n+2);
end


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


X0=[1, 1, 0, 0, 20, 1, 0, 0, 1, 0];
P0=eye(n)*2; P0(3,3)=10; P0(4,4)=10; P0(n, n)=0;

nZk=size(Zk,2);

% --- filter f_update (measurement model)
H=zeros(m,n); 

H(1,1)=1/calib.pix2cm(1); H(2,2)=1/calib.pix2cm(2); H(1:2,6)=calib.center;
H([3 4 5 6 7],[5 7 8 9 10])=eye(5); Hk=H;
R=eye(m); R(1:2,1:2)=diag(1./calib.pix2cm)*2; % note how we put the error in terms of the number of pixels per cm
R(3,3)=20; R(4:7,4:7)=0; % measurement noise sqrt of these values is measurement noise
Hinv=@(z) [(z(1:2)-calib.center).*(ones(2,1).*calib.pix2cm'); 0; 0; z(3); 1; z(4:m)]; % converts from Z to X

f_update=@(k, as1, Xh, P,  Zk) kalmanUpdateMT(k, as1, Xh, P, Zk(1:m,:), Hk, R, n);
initialize=@(k, Zk, idnew) init(k, Zk(1:m,:), Hinv, X0, P0, n, ci, idnew);

% --- f_predict (motion model)
[F, Fk]=motion_model(dt, n);
Q=diag([2 2 10 10 5 0 1 1 1 1]); %sqrt of these values is the std for the process noise; 0 in the last column implies that this value stays constant
f_predict=@(k, Xh, P) kalmanPredictMT(k, Xh, P, F, Q, n);

% --- associate
switch da_type
    case {0, 'gnn'}
        da1=@(D) gnnda(D, gt); 
    case {1, 'gnn+nnda'}
        da1=@(D) gnnda_nnda(D, gt); 
    case {3, 'nnda'}
        da1=@(D) nnda(D);
end

% --- measurement model
zhat=@(k, r, Xh) Hk*Xh(r, k);
zh_cov=@(c,r,P) Hk*P(r,c)*Hk' + R;

% --- measurement covariance
gvfun=@(Zh, Zk, S) (Zh-Zk)'/S*(Zh-Zk);
costmatrix=@(Zh, Zk, S) costmatrix2d(Zh, Zk, S, gvfun);

Zh=zeros(m,nt);
S=zeros(m,m,nt);

% prediction / motion model
if k>1
    [Xh, P]=f_predict(k, Xh, P);
end
for ii=1:nt
    [r c]=getind(n, 1, ii, 1:n, n);
    if Xh(r(1),1)
        Zh(:,ii)=zhat(1, r, Xh);
        S(:,:,ii)=zh_cov(c, r, P);
    end
end
Zh=Zh(:,Zh(1,:)~=0);

% matching measurements to targets
if ~isempty(Zk)
    D = costmatrix(Zh, Zk, S);

    % associate 
    as1=da1(D);

    % f_update
    [Xh P]=f_update(k, as1, Xh, P, Zk);
else % if there is no measurement then terminate the target
    as1=0;
    Zk=zeros(m,1);
    % f_update
    [Xh P]=f_update(k, as1, Xh, P, Zk);
end

% initialize
if ~isempty(Zk)
	Zk=Zk(:,~ismember(1:nZk,as1));
end
[Xi Pi]=initialize(k, Zk, idnew);


% put it back in the output form of X
Xru(:,ci+1:n+ci)=reshape(Xh, n, sum(curr))';
Xru(:,ck)=k;
Xru(:,cid)=Xr(curr,cid);

% put it back in the output form for covariance
for ii=1:sum(curr)
    P1=P(n*(ii-1)+1:n*ii,1:n);
    Prt(ii,:)=P1(:);
end
Pru=[];
if sum(curr)
    Pru(:,ci+1:n*n+ci)=Prt;
    Pru(:,ck)=k;
    Pru(:,cid)=Pr(curr,cid);
end


Xr=[Xr; Xru; Xi];
Pr=[Pr; Pru; Pi];

% delete all records that are trklen steps ago, current step, future, and
% not updated
delk=Xr(:,ck)<=k-trklen | Xr(:,ck)>k | Xr(:,3)==0;
Xr(delk,:)=[];
Pr(delk,:)=[];


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

function [Xh P] = kalmanUpdate(Xh, P, Z, H, R)
% function [Xh P] = kalmanUpdate(Xh, P, Z, Klmn)
%
% H is H
% R is R

% get state vector
nx = size(Xh,1);

% measurement covariance
S = H*P*H' + R;

% projected estimate
Zh = H*Xh;

% innovation
inov = Z-Zh;

% Klmn gain
W = P*H'/S;

[msg msgid]=lastwarn;
if strcmp(msgid, 'MATLAB:illConditionedMatrix') || strcmp(msgid, 'MATLAB:singularMatrix')
%     keyboard;
end 

% state update
Xh = Xh + W*inov;

% covariance update
P = (eye(nx)-W*H)*P;


function [Xh P]=kalmanUpdateMT(k, as1, Xh, P, Zk, Hk, R, n)

k=1;

% update
for tt=1:size(Xh,1)/n
    [r c]=getind(n, k, tt, 1:n, n);
    if numel(as1) >= tt
        if as1(tt)
            [Xh(r,k), P(r,c)]=kalmanUpdate(Xh(r,k), P(r,c), Zk(:,as1(tt)),  Hk, R);
        else
            Xh(r,k)=0;
            P(r,c)=0;
        end
    end
end

function [Xh, P] = kalmanPredict(Xh, P, F, Q)
%function [Xh, P] = kalmanPredict(Xh, P, Klmn)
%
% F is nxn transition matrix
% Q is nxn disturbance matrix


% state
Xh= F*Xh;

% covariance
P = F*P*F' + Q;

function [Xh, P]= kalmanPredictMT(k, Xh, P, F, Q, n)
k=1;
for ii=1:size(Xh,1)/n
    [r c]=getind(n, k, ii, 1:n, n);
    if Xh(r(1),1)
        [Xh(r,k), P(r,c)]=kalmanPredict(Xh(r,1), P(r,c), F, Q);
    end
end

function [Xi Pi]=init(k, Zk, Hinv, X0, P0, n, ci, idnew)
nZk=size(Zk,2);
Xi=zeros(nZk,n+ci);
Pi=zeros(nZk,n*n+ci);

% initialize
if ~isempty(Zk)
    for zz=1:nZk
        Xi(zz,:)=[k, idnew,  X0]; 
        Pi(zz,:)=[k, idnew, P0(:)'];
        Xi(zz, ci+1:n+ci)=Hinv(Zk(:,zz));
        idnew=idnew+1;
    end
end


function as1 = nnda(D)

[jnk as1]=min(D,[],1);

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


