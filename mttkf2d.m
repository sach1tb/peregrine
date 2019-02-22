function [Xr, Pr]=mttkf2d(Xr, Pr,  Zk, dt, calib, target, gatingThreshold, mode)
%function [Xr Pr]=mttkf2d(Xr, Pr,  Zk, dt, calib, target, mode)
%
% 2D kalman filter
% Xr is the [* x 6+2] row matrix for the current time step from the csv file
% Pr is the [* x 36+2] corresponding rows that will make the covariance matrix
% Zk is the [3 x *] matrix measurement for the current time step
% dt is the time between frames (e.g. 1/30 seconds)
% calib is the calibration matrix for the camera
% mode is for data association, set it to 1 for repair

global k
global trklen

trklen=5;

nZ=3; 
nX=6;
ck=1;
cid=2;
ci=2;

%%%%
if isempty(Xr)
    Xr=zeros(1,nX+ci);
    Pr=zeros(1,nX^2+ci);
end

%%% for adding a point
if size(Xr,1)==1 && Xr(1)
    P0=eye(nX)*2; P0(3,3)=10; P0(4,4)=10; P0(nX, nX)=0;
    Pr=[k Xr(2), P0(:)'];
end

if size(Xr,2) > nX+2
    Xr=Xr(:,1:nX+2);
end


curr=(Xr(:,ck)==k-1 & Xr(:,ck)~=0);

Xk=Xr(curr,ci+1:nX+2)'; % this is just the current time-step
Xh=Xk(:);
Pk=Pr(curr,ci+1:nX^2+2);
P=[];
for ii=1:sum(curr)
    P(nX*(ii-1)+1:nX*ii,1:nX)=reshape(Pk(ii,:),nX,nX);
end

nt=sum(curr);
ids=unique(Xr(:,cid));
if isempty(ids)
    idnew=1;
else
    idnew=max(ids)+1;
end

switch target
    case 2 % zebrafish
        % variance cm/s 
        Q=diag([0 0 9 9 5 0]); 
        X0=[1, 1, 0, 0, 20, 1];
        % average zebrafish speed is 5 cm/s
        P0=eye(nX)*0; P0(3,3)=25; P0(4,4)=25;
   case 3 % mosquitos
        % variance cm/s 
        Q=diag([0 0 5 5 5 0]); 
        X0=[1, 1, 0, 0, 20, 1];
        % average zebrafish speed is 3 cm/s
        P0=eye(nX)*0; P0(3,3)=9; P0(4,4)=9;    
    case 4 % pillbugs
        Q=diag([0 0 1 1 5 0]); % variance cm/s 
        X0=[1, 1, 0, 0, 20, 1];
        P0=eye(nX)*0; P0(3,3)=1; P0(4,4)=1; 
    case 1 % danios
        Q=diag([0 0 9 9 5 0]); % variance cm/s 
        X0=[1, 1, 0, 0, 20, 1];
        % average danio speed is 6 cm/s
        P0=eye(nX)*0; P0(3,3)=36; P0(4,4)=36; P0(nX, nX)=0;
end

nZk=size(Zk,2);

% --- filter f_update (measurement model)
H=zeros(nZ,nX); 

H(1,1)=1/calib.pix2cm(1); H(2,2)=1/calib.pix2cm(2); H(1:2,nX)=calib.center;
H(3,5)=1; Hk=H;
% note how we put the error in terms of the number of pixels per cm
R=eye(nZ); R(1:2,1:2)=diag(1./calib.pix2cm)*2; 
% measurement noise sqrt of these values is measurement noise
R(nZ,nZ)=20;
% converts from Z to X
Hinv=@(z) [(z(1:2)-calib.center).*(ones(2,1).*calib.pix2cm'); 0; 0; z(3); 1]; 

f_update=@(k, as1, Xh, P,  Zk) kalmanUpdateMT(k, as1, Xh, P, Zk(1:nZ,:), ...
    Hk, R, nX);
initialize=@(k, Zk, idnew) init(k, Zk(1:nZ,:), Hinv, X0, P0, nX, ci, idnew);

% --- f_predict (motion model)
[F, ~]=motion_model(dt, nX);
% sqrt of these values is the std for the process noise; 
% 0 in the last column implies that this value stays constant
f_predict=@(k, Xh, P) kalmanPredictMT(k, Xh, P, F, Q, nX);

% --- associate
switch mode
    case {1, 'gnn'}
        da1=@(D) gnnda(D, gatingThreshold);  
    case {2, 'nnda'}
        da1=@(D) nnda(D);
end

% --- measurement model
zhat=@(k, r, Xh) Hk*Xh(r, k);
zh_cov=@(c,r,P) Hk*P(r,c)*Hk' + R;

% --- measurement covariance
gvfun=@(Zh, Zk, S) (Zh-Zk)'/S*(Zh-Zk);
costmatrix=@(Zh, Zk, S) costmatrix2d(Zh, Zk, S, gvfun);

Zh=zeros(nZ,nt);
S=zeros(nZ,nZ,nt);

% prediction / motion model
if k>1
    [Xh, P]=f_predict(k, Xh, P);
end
for ii=1:nt
    [r, c]=getind(nX, 1, ii, 1:nX, nX);
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
    [Xh, P]=f_update(k, as1, Xh, P, Zk);
else % if there is no measurement then terminate the target
    as1=0;
    Zk=zeros(nZ,1);
    % f_update
    [Xh, P]=f_update(k, as1, Xh, P, Zk);
end

% initialize
if ~isempty(Zk)
	Zk=Zk(:,~ismember(1:nZk,as1));
end
[Xi, Pi]=initialize(k, Zk, idnew);


% put it back in the output form of X
Xru(:,ci+1:nX+ci)=reshape(Xh, nX, sum(curr))';
Xru(:,ck)=k;
Xru(:,cid)=Xr(curr,cid);

% put it back in the output form for covariance
for ii=1:sum(curr)
    P1=P(nX*(ii-1)+1:nX*ii,1:nX);
    Prt(ii,:)=P1(:);
end
Pru=[];
if sum(curr)
    Pru(:,ci+1:nX*nX+ci)=Prt;
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
        % position only
        D(zz,tt)=gvfun(Zh(1:2,tt), Zk(1:2,zz), S(1:2,1:2,tt)); 
        % position and size
%         D(zz,tt)=gvfun(Zh(:,tt), Zk(:,zz), S(:,:,tt)); 
    end
end

function [F, Fk]= motion_model(dt,nX)
Fk=eye(nX); Fk(1,3)=dt; Fk(2,4)=dt; 
F=Fk; 

function [Xh, P] = kalmanUpdate(Xh, P, Z, H, R)
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

[~, msgid]=lastwarn;
if strcmp(msgid, 'MATLAB:illConditionedMatrix') || ...
        strcmp(msgid, 'MATLAB:singularMatrix')
%     keyboard;
end 

% state update
Xh = Xh + W*inov;

% covariance update
P = (eye(nx)-W*H)*P;


function [Xh, P]=kalmanUpdateMT(~, as1, Xh, P, Zk, Hk, R, nX)

k=1;

% update
for tt=1:size(Xh,1)/nX
    [r, c]=getind(nX, k, tt, 1:nX, nX);
    if numel(as1) >= tt
        if as1(tt)
            [Xh(r,k), P(r,c)]=kalmanUpdate(Xh(r,k), P(r,c), ...
                Zk(:,as1(tt)),  Hk, R);
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

function [Xh, P]= kalmanPredictMT(~, Xh, P, F, Q, nX)

k=1;
for ii=1:size(Xh,1)/nX
    [r, c]=getind(nX, k, ii, 1:nX, nX);
    if Xh(r(1),1)
        [Xh(r,k), P(r,c)]=kalmanPredict(Xh(r,1), P(r,c), F, Q);
    end
end

function [Xi, Pi]=init(k, Zk, Hinv, X0, P0, nX, ci, idnew)
nZk=size(Zk,2);
Xi=zeros(nZk,nX+ci);
Pi=zeros(nZk,nX*nX+ci);

% initialize
if ~isempty(Zk)
    for zz=1:nZk
        Xi(zz,:)=[k, idnew,  X0]; 
        Pi(zz,:)=[k, idnew, P0(:)'];
        Xi(zz, ci+1:nX+ci)=Hinv(Zk(:,zz));
        idnew=idnew+1;
    end
end


function as1 = nnda(D)

[~, as1]=min(D,[],1);

function as1 = gnnda(D, thresh)

as=munkres(D);
idx=find(as>0);

% remove all associations beyond a certain threshold
% as is number of targets with each having a measurement number
as(idx(D(idx)>thresh))=0; 

[val, as1]=max(as,[],1);
as1(val==0)=0;



function as1 = gnnda_nnda(D, thresh)


as=munkres(D);
idx=find(as>0);

% remove all associations beyond a certain threshold
% as is number of targets with each having a measurement number
as(idx(D(idx)>thresh))=0; 

[val, as1]=max(as,[],1);
as1(val==0)=0;

%%%% nearest neighbor
for ii=find(as1==0)
    [vv, zz]=min(D(:,ii));
    if vv<thresh/2
        as1(ii)=zz;
    end
end


