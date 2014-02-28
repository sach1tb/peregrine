function get_obs(type, varargin)
%function get_obs(type, 'Parameter Name', Parameter value, ...)
%
% type is the observable:
% 'Nearest neighbor distance'
% 'Polarization'
% 'Speed'
% 'Turn rate'
% 'Freezing'
% 'Thrashing'
% 'Swimming'
% 'Excursions'
% 'Time spent left'
% 'Time spent right'
% 'Acceleration'
% 'Distance from stimulus'
% 'Time far from stimulus'
% 'Time near stimulus'
% 'Number of entries into a region'
% 'Leadership time lag xcorr'
% 'Transfer entropy 1-2'
%
% Other parameters that are needed based on the observable
%
% 'Fps' = frame rate
% 'FocalId' = id of the fish if you want to analyse a single fish. if not
% specified, group average is reported
% 'NumTargets' = number of targets
% 'TrialTimeMin' = trial time in minutes. if this is a 1x2 vector then only
% the time between those two values in minutes is extracted
% 'RoiCm' = roi is a 1x2 vector e.g. [56 30]. a 1x4 vector assumes exact
% roi, e.g. [-27 -15 27 15] (xmin ymin xmax ymax]
% 'StimPosCm'= position of stimulus if only one value is specified then
% distance along the same is computed e.g. [10] or [5 10] 
% 'RegionBounds' = Boundaries of a region to compute number of entries;
% specify min and maximum e.g. [15, 28], or [-28, -15]
% 'XcorrIds' = ids of the targets between whom you wish to find
% xcorrelaion e.g [1 2], if the second id is 0 the centroid of all other
% ids is used e.g. [1 0]
% 'TrajSec' = segments of trajectory in seconds that must be used to x corr
% 'TankPartitions' = number of partitions you want to divide the tank into
%
% Examples:
% 1) find the percentage time spent freezing by fish identified by id 1
% get_obs('Freezing', 'fps', 24, 'numTargets', 5,'TrialTimeMin', 5, 'Focalid', 1)
%
% 2) find the percentage time spent freezing by all fish
% get_obs('Freezing', 'fps', 24, 'numTargets', 5,'TrialTimeMin', 5)
%
% 3) find the average speed of all fish
% get_obs('Speed', 'fps', 24, 'numTargets', 5,'TrialTimeMin', 5)
%
% 4) find the distance to a stimulus placed at the right side of a [70x30]
% cm region [CARE should be taken that the roi does not include stimulus
% region] NOTE: the position of stimulus is half that of the RoiCm(1) on
% each side
% get_obs('Distance from stimulus', 'fps', 24, 'numTargets', 5, 'TrialTimeMin', 5, 'RoiCm', [70 30], 'StimPosCm', 35)
% for left side
% get_obs('Distance from stimulus', 'fps', 24, 'numTargets', 5, 'TrialTimeMin', 5, 'RoiCm', [70 30], 'StimPosCm', -35)
%
% 5) find the number of entries into a region
% get_obs('Number of entries into a region', 'fps', 15, 'numTargets',1, 'TrialTimeMin', 5, 'FocalId', 1, 'RegionBounds', [15, 28])
addpath('../');

% parse inputs
p=inputParser;

p.addParamValue('Fps', 0);
p.addParamValue('FocalId', 0);
p.addParamValue('NumTargets', nan);
p.addParamValue('TrialTimeMin', 0);
p.addParamValue('RoiCm', [0 0]);
p.addParamValue('StimPosCm', [nan]);
p.addParamValue('RegionBounds', [0 0]);
p.addParamValue('XcorrIds', [0 0]);
p.addParamValue('TrajSec', 0);
p.addParamValue('TankPartitions', 3);

p.parse(varargin{:});

fid=p.Results.FocalId;
fps=p.Results.Fps;
nt=p.Results.NumTargets;
nmin=p.Results.TrialTimeMin;
roi=p.Results.RoiCm;
stimpos=p.Results.StimPosCm;
regionbounds=p.Results.RegionBounds;
xcorrids=p.Results.XcorrIds;
trajsec=p.Results.TrajSec;
tparts=p.Results.TankPartitions;

%%%% processing
if fid && nt ~=1
    warning('MATLAB:numtargets','Number of expected targets should be 1 for a focal fish');
end

% validating
if ~fps
    error('Fps must be specified');
end

if ~nmin
    error('TrialTimeMin must be specified');
end

if isnan(stimpos) && strcmp('Distance from stimulus', type)
    error('StimPosCm must be specified for Distance from stimulus')
end

if isnan(nt) && strcmp('Nearest Neighbor distance', type)
    error('NumTargets must be specified for Nearest Neighbor Distance')
end

if numel(fid)~=2 && strcmp('Transfer entropy 1-2', type)
    error('Specify two focal ids with transfer entropy from 1 to 2')
end

if (isnan(stimpos) && strcmp('Time near stimulus', type)) || ...
   (isnan(stimpos) && strcmp('Time far from stimulus', type))
        error('StimPosCm must be specified for Time near/far from stimulus')
end

if ~sum(regionbounds) && strcmp('Number of entries into a region', type)
    error('Specify the region boundaries, min and maximum [15, 28], or [-28, -15]');
end

fpm=fps*60; % frames per minute = frames per second * 60 seconds
nfrm=fpm*nmin;

[FileName,PathName]=uigetfile('*.csv', 'Select the data file(s)', 'MultiSelect', 'on');

if ~iscell(FileName)
    FileName=cellstr(FileName);
end

switch type
    case 'Nearest neighbor distance'
        bin_centers=0:1.25:25; units='(cm)'; acr='annd'; 
        get_data=@(PathName, FileName) get_annd(PathName, FileName, nt, nfrm);
    case 'Polarization'
        bin_centers=0:0.05:1; units=''; acr='pol'; 
        get_data=@(PathName, FileName) get_pol(PathName, FileName, nt, nfrm);
    case 'Speed'
        bin_centers=0:1:10; units='(cm/s)'; acr='speed'; 
        get_data=@(PathName, FileName) get_speed(PathName, FileName, nt, nfrm, fid);
    case 'Turn rate'
        bin_centers=0:5:200; units='(degree/s)'; acr='turnrate'; 
        get_data=@(PathName, FileName) get_turn_rate(PathName, FileName, nt, nfrm, fid, fps);        
    case 'Freezing'
        bin_centers=0:.1:1; units='(%)'; acr='ptf'; acc_t=120; wall_delta=3; tsec=2; rad=2;
        get_data=@(PathName, FileName) get_freeze(PathName, FileName, nt, nfrm, fid, fps,...
                                                    tsec, acc_t, roi, wall_delta, rad); 
    case 'Thrashing'
        bin_centers=0:.1:1; units='(%)'; acr='ptt'; acc_t=120; wall_delta=3; tsec=2; rad=2;
        get_data=@(PathName, FileName) get_thrash(PathName, FileName, nt, nfrm, fid, fps,...
                                                    tsec, acc_t, roi, wall_delta, rad);
    case 'Swimming'
        bin_centers=0:.1:1; units='(%)'; acr='ptt'; acc_t=120; wall_delta=3; tsec=2; rad=2;
        get_data=@(PathName, FileName) get_swimm(PathName, FileName, nt, nfrm, fid, fps,...
                                                    tsec, acc_t, roi, wall_delta, rad);                                             
    case 'Excursions'
        bin_centers=0:.1:1; units='(%)'; acr='excur'; 
        get_data=@(PathName, FileName) get_excursions(PathName, FileName, nt, nfrm, fid);    
    case 'Leadership time lag xcorr'
        bin_centers=-.5:.1:.5; units='(sec)'; acr='timelag';
        get_data=@(Pathname, FileName) get_leadership_xcorr(PathName, FileName, nt, nfrm, fps, trajsec, xcorrids);
    case 'Time spent left'
        bin_centers=0:.1:1; units='(%)'; acr='tsl'; 
        get_data=@(PathName, FileName) get_tsl(PathName, FileName, nt, nfrm, fid); 
    case 'Time spent right'
        bin_centers=0:.1:1; units='(%)'; acr='tsr'; 
        get_data=@(PathName, FileName) get_tsr(PathName, FileName, nt, nfrm, fid);    
    case 'Acceleration'
        bin_centers=0:1:20; units='(cm/s^2)'; acr='acc'; 
        get_data=@(PathName, FileName) get_accel(PathName, FileName,nt, nfrm, fid, fps);
    case 'Distance from stimulus'
        bin_centers=0:2:100; units='(cm)'; acr='d2s'; 
        get_data=@(PathName, FileName) get_dist2stim(PathName, FileName, nt, nfrm, fid, stimpos);  
    case 'Time near stimulus'
        bin_centers=0:.1:1; units='(%)'; acr='tsns'; 
        get_data=@(PathName, FileName) get_time_near_stim(PathName, FileName, nt, nfrm, fid, stimpos, roi, tparts);  
    case 'Time far from stimulus'
        bin_centers=0:.1:1; units='(%)'; acr='tsfs'; 
        get_data=@(PathName, FileName) get_time_far_stim(PathName, FileName, nt, nfrm, fid, stimpos, roi, tparts);
    case 'Number of entries into a region'
        bin_centers=0:2:50; units=''; acr='freqr'; 
        get_data=@(PathName, FileName) get_freq_into_region(PathName, FileName, nt, nfrm, fid, regionbounds);
    case 'Transfer entropy 1-2'
        bin_centers=0:.2:2; units=''; acr='te12'; 
        get_data=@(PathName, FileName) get_transfer_entropy_12(PathName, FileName, nt, nfrm, fid);
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% anal='Variation in '; 
% anal='Average ';

% figure(1); gcf; clf;
% ha=tight_subplot(1, numel(cselect), [0.025, 0.01], [0.2, 0.025], [0.05 .01]);


% ad_test_1way=nan(10, max(cselect));
% ad_test_2way=nan(10, 5, max(cselect));

[data permin]=get_data(PathName, FileName);
% cut the last five minutes
fprintf('\nAverage *%s* %s \n', type, units);
stats=nanmean(data,2);
fprintf('FileName\t\t Avg. Value \t\t P1 \t P2 \t P3 \t P4 \t P5\n', type, units);
fprintf('------------------------------------------------\n');
for jj=1:size(FileName,2)
    fprintf('%s\t\t %.3f \t\t %.3f \t %.3f \t %.3f \t %.3f \t %.3f\n', FileName{1,jj}, stats(jj), permin(jj,1), permin(jj,2), ...
                                                                        permin(jj,3), permin(jj,4), permin(jj,5));
end

uav=~isnan(data);
fprintf('%% data used = %.2f +/- %.2f\n', mean(mean(uav,2))*100, std(mean(uav,2)*100));

% stats
resp=input('Do you want to run One-way ANOVA stats on this data? []=yes, other=no: ');

if isempty(resp)
    for jj=1:size(FileName,2)
        tmp=textscan(FileName{jj}, '%s', 'Delimiter', '_');
        cond{jj}=tmp{1}{2};
    end

    fprintf('The stats are run with the following conditions as per the naming of files: \n')
    cc=unique(cond);
    for ii=1:size(cc,2)
        fprintf('%s\n',cc{ii});
    end
    figure(1); gcf; clf;
    subplot(1,2,1); gca;
    mystat('anovan', stats, cond, '', 3);
    
    subplot(1,2,2); gca;
    plot(permin');
end

function flag=filter_nt(rf, nt)
flag=double(size(rf,2)==nt);

% Nearest Neighbor Distance
function [dat permin]=get_annd(PathName, FileName, nt, nfrm)
id=0;
if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end
for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    for k=1:size(data.X,2)
        rf =processk(data.X, Xi, k, id);
        if filter_nt(rf,nt)
            dat(ii,k)=annd(rf(1:2,:));
        else
            dat(ii,k)=nan;
        end
    end
end

permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
if size(permin,2)==1, permin=permin'; end

function [dat perm]=annd(r)
% function annd(r)
% average nearest-neighbor distance
%
% r is a d x n matrix with d dimensions and n targets

[d n]=size(r);
D=eye(n)*100000000;

for ii=1:n
    for jj=ii+1:n
        D(ii,jj)=norm(r(:,ii)-r(:,jj));
    end
end

D1=D+D';
D=D1;

[perm idx]=min(D);

dat=mean(perm);


% Speed
function [dat permin]=get_speed(PathName, FileName, nt, nfrm, id)

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);

     if id
        r=getind(Xi.nX, 1, id, 1:Xi.nX,1);
        data.X=data.X(r,:);
        nt=1;
     end
    

    v1=data.X(3:Xi.nX:end,:);
    v2=data.X(4:Xi.nX:end,:);
   
    nz=sum(v1~=0,1);
   
    dat(ii,1:size(v1,2))=sum(sqrt(v1.^2+v2.^2),1)./nz;
    dat(ii,nz~=nt)=nan;
end

permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
if size(permin,2)==1, permin=permin'; end

% Turn rate
function [dat permin]=get_turn_rate(PathName, FileName, nt, nfrm, id, fps)

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end
for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);

     if id
        r=getind(Xi.nX, 1, id, 1:Xi.nX,1);
        data.X=data.X(r,:);
        nt=1;
     end
    
    
    v1=data.X(3:Xi.nX:end,:);
    v2=data.X(4:Xi.nX:end,:);
   
    % get current timestep and previous timestep
    v1k=v1(:,1:end-1);
    v1km1=v1(:,2:end);
    
    v2k=v2(:,1:end-1);
    v2km1=v2(:,2:end);
    
    % normalize
    v1kn=v1k./sqrt(v1k.^2+v2k.^2);
    v2kn=v2k./sqrt(v1k.^2+v2k.^2);
    
    v1km1n=v1km1./sqrt(v1km1.^2+v2km1.^2);
    v2km1n=v2km1./sqrt(v1km1.^2+v2km1.^2);
    
    dp=v1kn.*v1km1n+v2kn.*v2km1n;
    
    nz=sum(v1k~=0,1);
   
    tr=acos(dp)*fps*180/pi;
    tr(v1k==0 & v2k==0)=0;
    dat(ii,1:size(tr,2))=sum(real(tr),1)./nz;
    dat(ii,nz~=nt)=nan;
end

permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
if size(permin,2)==1, permin=permin'; end

% Polarization
function [dat permin]=get_pol(PathName, FileName, nt, nfrm)

id=0;
if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end
for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    for k=1:size(data.X,2)
        [rf vf]=processk(data.X, Xi, k, id);
        if filter_nt(rf,nt)
            dat(ii,k)=pol(vf);
        else
            dat(ii,k)=nan; 
        end
    end
end

permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
if size(permin,2)==1, permin=permin'; end

% Freezing 
function [dat permin]=get_freeze(PathName, FileName, nt, nfrm, id, fps, tsec, acc_t, roi, wall_delta, rad)

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    if id
        data.X=data.X(Xi.nX*(id-1)+1:Xi.nX*id,:);
    end
    
    % frames corresponding to 2 sec and 2 cm radius
    fr1=locomotory_behavior(data.X, Xi.nX, fps, tsec, acc_t, roi, wall_delta, rad, Xi); 
    
    dat(ii,1:size(fr1,2))=fr1;
end

permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
if size(permin,2)==1, permin=permin'; end

function [freeze thrash swimm]=locomotory_behavior(X, nX, fps, tsec, acc_t, roi, wall_delta, dist_t, Xi)
% function freeze=locomotory_behavior(X,nX,n,dist_t)
%
% fr1=locomotory_behavior(data.X, Xi.nX, 48, 2); % frames corresponding to 2 sec and 2 cm radius
%
% This script will mark all frames that are frozen as 1
%

n=tsec*fps;

tf=size(X,2);
freeze=zeros(1,tf);
thrash=freeze;


r1=X(1:nX:end,:);
r2=X(2:nX:end,:);


 %%% remove thrashing from these frames
% velocities
v1=X(3:Xi.nX:end,:);
v2=X(4:Xi.nX:end,:);

% accelerations
a1=diff(v1, [], 2)*fps;
a2=diff(v2, [], 2)*fps;

% magnitude
accel=sqrt(a1.^2 + a2.^2);


%%%%%% [1] V. Kopman, J. Laut, G. Polverino, and M. Porfiri, 
% "Closed-loop control of zebrafish response using a bioinspired robotic-fish in a preference test",
% Journal of the Royal Society Interface, vol. 10, no. 78, p. 20120540, 2013.

for ff=1:size(r1,1)
    for ii=1:n:size(r1,2)-n
        rt=[r1(ff,ii:ii+n-1);
            r2(ff,ii:ii+n-1);];
        dist=zeros(1, size(rt,2)-1);
        for jj=2:size(rt,2)
            dist(1,jj-1)=norm(rt(:,1)- rt(:,jj));
        end
        if max(dist(1,:))<dist_t && ~sum(dist==0)
            if numel(roi)==2
                if sum(accel(ff,ii:ii+n-1)) > acc_t && ... % if the fish is accelerating above a value 
                    (sum(abs(rt(1,:))>roi(1)/2-wall_delta) || ... % if the fish is near the wall
                    sum(abs(rt(2,:))>roi(2)/2-wall_delta))
                    thrash(ii:ii+n-1)=1;
                else
                    freeze(ii:ii+n-1)=1;
                end
            elseif numel(roi)==4
                if sum(accel(ff,ii:ii+n-1)) > acc_t && ...
                    (sum(rt(1,:)<roi(1)+wall_delta) || ...
                     sum(rt(2,:)<roi(2)+wall_delta) || ...
                     sum(rt(1,:)>roi(3)-wall_delta) || ...
                     sum(rt(2,:)>roi(4)-wall_delta))
                    thrash(ii:ii+n-1)=1;
                else
                    freeze(ii:ii+n-1)=1;
                end
            end
        end
    end
end

swimm=1-(freeze+thrash);

%%% faster version
% for ff=1:size(r1,1) % number of fish
%         r1t=zeros(n-1,tf-n+1);
%         r2t=r1t;
%         for jj=1:n-1
%             r1t(jj,:)=r1(ff,jj+1:tf-n+jj+1)-r1(ff,1:tf-n+1);
%             r2t(jj,:)=r2(ff,jj+1:tf-n+jj+1)-r2(ff,1:tf-n+1);
%         end
%         rchk=sqrt(r1t.^2+r2t.^2);
%         rchk=max(rchk,[],1);
%         
%         freeze(rchk<dist_t & rchk>0)=1;  
% end


% even faster....
%
% nf=size(r1,1);
% r1t=zeros(nf, tf-n+1, n-1);
% r2t=r1t;
% for jj=1:n-1
%     r1t(:,:,jj)=r1(:,jj+1:tf-n+jj+1)-r1(:,1:tf-n+1);
%     r2t(:,:,jj)=r2(:,jj+1:tf-n+jj+1)-r2(:,1:tf-n+1);
% end
% rchk=sqrt(r1t.^2+r2t.^2);
% rchk=max(rchk, [], 3); % maximum distance in the time
% rchk(rchk==0)=dist_t; % to remove zero entries
% rchk=min(rchk, [], 1); % any fish
% 
% freeze(rchk<dist_t & rchk>0)=1; 
% toc


% Thrashing 
function [dat permin]=get_thrash(PathName, FileName, nt, nfrm, id, fps, tsec, acc_t, roi, wall_delta, rad)

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    if id
        data.X=data.X(Xi.nX*(id-1)+1:Xi.nX*id,:);
    end
    
    % frames corresponding to 2 sec and 2 cm radius
    [fr1 tr1]=locomotory_behavior(data.X, Xi.nX, fps, tsec, acc_t, roi, wall_delta, rad, Xi); 

    dat(ii,1:size(tr1,2))=tr1;
end

permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
if size(permin,2)==1, permin=permin'; end

% Swimming 
function dat=get_swimm(PathName, FileName, nt, nfrm, id, fps, tsec, acc_t, roi, wall_delta, rad)

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    if id
        data.X=data.X(Xi.nX*(id-1)+1:Xi.nX*id,:);
    end
    
    % frames corresponding to 2 sec and 2 cm radius
    [fr1 tr1 sw1]=locomotory_behavior(data.X, Xi.nX, fps, tsec, acc_t, roi, wall_delta, rad, Xi); 

    dat(ii,1:size(sw1,2))=sw1;
end

permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
if size(permin,2)==1, permin=permin'; end

% Excursions
function dat=get_excursions(PathName, FileName, nt, nfrm, id)

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    for k=1:size(data.X,2)
        rf=processk(data.X, Xi, k, 0);
        [jnk perm]=annd(rf);
        
        if filter_nt(rf,nt)
            annd_pfish(:,k)=perm';
        else
            annd_pfish(:,k)=zeros(nt,1);
        end
    end
    dat(ii,:)=shoal_membership_f1(annd_pfish, id);
end

function excurf1=shoal_membership_f1(annd_pfish, id)
%
% ensure that the first fish is which you want to compute membership of
%
% ref:
% [1] N. Miller and R. Gerlai, ?Redefining membership in animal groups.,? 
% Behavior Research Methods, vol. 43, no. 4, pp. 964?970, Dec. 2011.


% fig. 1 (a) 

ts_annd_f1=annd_pfish(id,:);

% fig. 1 (b)
[N bin]=hist(annd_pfish(annd_pfish~=0));

[val idx]=max(N);
modev=bin(idx);

exceed_mode=ts_annd_f1>modev;

j=0;
msv=[]; msk=[]; msmx=[];
for k=2:numel(ts_annd_f1)
    if exceed_mode(k)
        if ~exceed_mode(k-1) || k==2
            j=j+1;
        end
        msv=[msv ts_annd_f1(k)];
        msk=[msk k];
        msmx=max(msv);
        ms(j).v=msv;
        ms(j).k=msk;
        ms(j).mx=msmx;
    else
        if exceed_mode(k-1)
            msv=[]; msk=[]; msmx=[];
        end
    end
end

% fig. 1 (c)
ms1=cat(1,ms.mx);

p95=prctile(ms1, 95);
excurf1=zeros(1, size(annd_pfish,2));
for jj=1:size(ms,2)
    if ms(jj).mx > p95
        excurf1(ms(jj).k)=1;
    end
end

% see above on line 24
excurf1(1)=1;



% leadership using xcorr
function [dat permin]=get_leadership_xcorr(PathName, FileName, nt, nfrm, fps, trajsec, xcorrids)

if ~xcorrids(1)
    error('fish id for x correlation must be specified..');
end

if numel(nfrm)==1
    dat=nan(size(FileName,2), nfrm);
else
    dat=nan(size(FileName,2), nfrm(2)-nfrm(1));
end

trajsec=trajsec*fps;

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
    
    data.X= clean_and_smooth(data.X,Xi,nt);

    % for each section of the trajectory
    for k=1:trajsec:size(data.X,2)-trajsec
        % exposed fish
        re=getind(Xi.nX, k, xcorrids(1), 3:4,1);
        vve=normv(data.X(re,k:k+trajsec));
        dir_e=atan2(vve(2,:), vve(1,:));
        
        % rest of the shoal
        v1=data.X(3:Xi.nX:end,k:k+trajsec);
        v2=data.X(4:Xi.nX:end,k:k+trajsec);

        nz=sum(v1~=0,1);

        v1us=sum(v1(2:end,:),1)./(nz-1);
        v2us=sum(v2(2:end,:),1)./(nz-1);
        
        dir_us=atan2(v2us, v1us);
        
        [xx lags]=xcorr(dir_e, dir_us);
%         [xx lags]=xcorr(vve(1,:), v1us);
        [val idx]=max(xx);
        tau=lags(idx)/fps;
        
%         tau=lags(idx)/fps*val;
        dat(ii,k:k+trajsec)=tau;
    end
end

permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
if size(permin,2)==1, permin=permin'; end

% Time spent left
function [dat permin]=get_tsl(PathName, FileName, nt, nfrm, id)

if ~id
    error('fish id must be specified..');
end

if numel(nfrm)==1
    dat=nan(size(FileName,2), nfrm);
else
    dat=nan(size(FileName,2), nfrm(2)-nfrm(1));
end

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    data.X=data.X(Xi.nX*(id-1)+1:Xi.nX*id,:);
    
    r1=data.X(1,:);
    
    dat(ii,:)=double(r1 < 0 & r1~=0);
end

permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
if size(permin,2)==1, permin=permin'; end

% Time spent right
function [dat permin]=get_tsr(PathName, FileName, nt, nfrm, id)

if ~id
    error('fish id must be specified..');
end

if numel(nfrm)==1
    dat=nan(size(FileName,2), nfrm);
else
    dat=nan(size(FileName,2), nfrm(2)-nfrm(1));
end

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    data.X=data.X(Xi.nX*(id-1)+1:Xi.nX*id,:);
    
    r1=data.X(1,:);
    
    dat(ii,:)=double(r1 > 0 & r1~=0);
end

permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
if size(permin,2)==1, permin=permin'; end

% Acceleration
function [dat permin]=get_accel(PathName, FileName, nt, nfrm, id, fps)

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end
for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    if id
        data.X=data.X(Xi.nX*(id-1)+1:Xi.nX*id,:);
    end
    
    % velocities
    v1=data.X(3:Xi.nX:end,:);
    v2=data.X(4:Xi.nX:end,:);
    
    % accelerations
    a1=diff(v1, [], 2)*fps;
    a2=diff(v2, [], 2)*fps;
    
    % magnitude
    accel=sqrt(a1.^2 + a2.^2);
    
    % number of targets tracked in every frame
    nf=sum(data.X(1:Xi.nX:end,:)~=0,1);
    
    % mean acceleration
    dat(ii,1:size(accel,2))=sum(accel,1)./nf(1:size(accel,2));
    
    % remove all those that are not good
    dat(ii,nf~=nt)=nan;
end

permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
if size(permin,2)==1, permin=permin'; end

% Distance from stimulus
function [dat permin]=get_dist2stim(PathName, FileName, nt, nfrm, id, stimpos)

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    if id
        data.X=data.X(Xi.nX*(id-1)+1:Xi.nX*id,:);
    end
    
    % position
    r1=data.X(1:Xi.nX:end,:); r1(r1==0)=nan;
    r2=data.X(2:Xi.nX:end,:); r2(r2==0)=nan;
    
    if numel(stimpos)==1
        dat(ii,:)=nanmean(abs(stimpos(1)-r1),1);
    else
        dat(ii,:)=nanmean(sqrt((stimpos(1)-r1).^2+(stimpos(2)-r2).^2),1);
    end
end

permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
if size(permin,2)==1, permin=permin'; end

% Time near stimulus
function [dat permin]=get_time_near_stim(PathName, FileName, nt, nfrm, id, stimpos, roi, tparts)

if ~id
    error('fish id must be specified..');
end

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    data.X=data.X(Xi.nX*(id-1)+1:Xi.nX*id,:);
    
    r1=data.X(1,:);
    
    dat(ii,:)=double(abs(r1-stimpos(1)) < roi(1)/tparts);
end

permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
if size(permin,2)==1, permin=permin'; end

% Time far from stimulus
function [dat permin]=get_time_far_stim(PathName, FileName, nt, nfrm, id, stimpos, roi, tparts)

if ~id
    error('fish id must be specified..');
end

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    data.X=data.X(Xi.nX*(id-1)+1:Xi.nX*id,:);
    
    r1=data.X(1,:);
    
    dat(ii,:)=double(abs(r1-stimpos(1)) > (tparts-1)*roi(1)/tparts);
end

permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
if size(permin,2)==1, permin=permin'; end


function [dat permin]=get_freq_into_region(PathName, FileName, nt, nfrm, id, regionbounds)

if ~id
    error('fish id must be specified..');
end

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    data.X=data.X(Xi.nX*(id-1)+1:Xi.nX*id,:);
    
    r1=data.X(1,:);
    
    for k=2:numel(r1)
        if r1(k) > regionbounds(1) && r1(k-1) < regionbounds(1)
            dat(ii,k)=1;
        end
        if r1(k) < regionbounds(2) && r1(k-1) > regionbounds(2)
            dat(ii,k)=1;
        end
    end
end

permin=squeeze(sum(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
if size(permin,2)==1, permin=permin'; end

dat=sum(dat,2);


% function [dat permin]=get_transfer_entropy_12(PathName, FileName, nt, nfrm, id, regionbounds)
% 
% if numel(id)~=2
%     error('both fish ids must be specified..');
% end
% 
% if numel(nfrm)==1
%     dat=zeros(size(FileName,2), nfrm);
% else
%     dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
% end
% 
% for ii =1:size(FileName,2)
%     [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
%  
%     data.X= clean_and_smooth(data.X,Xi,nt);
%     
%     data.X=data.X(Xi.nX*(id-1)+1:Xi.nX*id,:);
%     
%     r1=data.X(1,:);
%     
%     for k=2:numel(r1)
%         if r1(k) > regionbounds(1) && r1(k-1) < regionbounds(1)
%             dat(ii,k)=1;
%         end
%         if r1(k) < regionbounds(2) && r1(k-1) > regionbounds(2)
%             dat(ii,k)=1;
%         end
%     end
% end
% 
% permin=squeeze(sum(reshape(dat, [size(dat,1), size(dat,2)/5, 5]), 2));
% if size(permin,2)==1, permin=permin'; end
% 
% dat=sum(dat,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%________________________________________________________________

function [Xh Xi]=reformat_csv(csvfile, nfrm)

Xhcsv=csvread(csvfile);

% clean the file to ensure that some noise is removed here only
max_t=max(Xhcsv(:,2));
if max_t > 100
    fprintf('[W] found more than 100 targets! Possibly noise\n');
    Xhcsv=Xhcsv(Xhcsv(:,2)<100,:);
    max_t=max(Xhcsv(:,2));
end

if sum(Xhcsv(:,10))
    Xi=strXi('shape');
else
    Xi=strXi;
end


max_k=max(Xhcsv(:,1));
Xh=zeros(Xi.nX*max_t, max_k);

if size(Xhcsv,1) < 1
    error('csv file is empty. Check...');
end

% format the Xh
% tic
% for ii=1:size(Xhcsv,1)
%     k=Xhcsv(ii,1);
%     id=Xhcsv(ii,2);
%     if id && k
%         [r c]=getind(Xi.nX, k, id, 1:Xi.nX, 1);
%         Xh(r,c)=Xhcsv(ii,3:3+Xi.nX-1)';
%     end
% end
% toc

% alternatel
% tic
ids=unique(Xhcsv(:,2)); ids=ids(ids~=0)';

for ii=ids
    r=getind(Xi.nX, 1, ii, 1:Xi.nX, 1);
    rows=(Xhcsv(:,2)==ii & Xhcsv(:,1)~=0);
    Xh_k=Xhcsv(rows,1);
    Xh(r,Xh_k)=Xhcsv(rows,3:3+Xi.nX-1)';
end
% toc

% shave off the last nfrms
ff=Xhcsv(Xhcsv(:,2)==1 & Xhcsv(:,1),1); % get all the frames for id==1
if numel(nfrm)==1
    if isempty(ff), ff=1; end
    ff=min(ff); % get the first frame
    nfr=nfrm;
else
    ff=nfrm(1);
    nfr=nfrm(2)-nfrm(1);
end
try
    if size(Xh,2) < ff+nfr
        fprintf('[W] Appending %d frames ...\n', ff+nfrm-size(Xh,2)+1);
        Xh=[Xh (Xh(:,1)*0)*ones(1,ff+nfr-size(Xh,2)+1)];
    end
    Xh=Xh(:,ff:ff+nfr-1);
catch ME
    fprintf('[!] %s (%d frames missing)\n', csvfile, (ff+nfr)-size(Xh,2));
end

function Xi=strXi(type)

if nargin ~=1
    type='pos';
end

switch type
    case 'pos'
                % state space 
        Xi.ri=1:2; % head position homogeneous coordinates
        Xi.rdi=3:4; % heading
        Xi.sz=5; % size or area
        %         Xi.or=6; % orientation
        Xi.nX=6;
    case 'shape'
                        % state space 
        Xi.ri=1:2; % head position homogeneous coordinates
        Xi.rdi=3:4; % heading
        Xi.sz=5; % size or area
        %         Xi.or=6; % orientation
        Xi.sh=7:8;
        Xi.hd=9:10;
        Xi.nX=10;
end

function X= clean_and_smooth(X,Xi, nt)

% remove all other tracks if nt tracks are more than 90% of the data
nt_k=sum(X(1:Xi.nX:end,:)~=0,2);
nfrm=size(X,2);
if sum(nt_k./nfrm > .9) == nt
    idx=find(nt_k./nfrm >.9);
    
    idx=((idx-1)*ones(1,Xi.nX)*Xi.nX+ones(numel(idx),1)*(1:Xi.nX))';
    idx=idx(:);
    
    X=X(idx,:);
end

win=5;
% 
for jj=1:Xi.nX:size(X,1)
    if sum(X(jj,:)~=0)<win
        X(jj:jj+Xi.nX-1,:)=0;
    end
end

for jj=1:size(X,1)
    X(jj,:)=sma(X(jj,:), win);
end


function [rf vf]=processk(X, Xi, k, id)

if id
    data=X(Xi.nX*(id-1)+1:Xi.nX*id,k);
    rf=data(1:2);
    if Xi.nX==10
        vf=data(9:10);
    else
        vf=data(3:4);
    end
else
    rf=[X(1:Xi.nX:end,k)'; X(2:Xi.nX:end,k)'];
    if Xi.nX==10
        vf=[X(9:Xi.nX:end,k)'; X(10:Xi.nX:end,k)'];
    else
        vf=[X(3:Xi.nX:end,k)'; X(4:Xi.nX:end,k)'];
    end
end
nz=rf(1,:)~=0;
rf=rf(:,nz) ;


vf=vf(:,nz) ; % no robot


if 0 % debug
    figure(1); gcf; clf;
    quiver(rf(1,:),rf(2,:), vf(1,:), vf(2,:), 'k');
    axis([-60 60 -60 60]);
end

