function get_obs(type, varargin)
%function get_obs(type, 'Parameter Name', Parameter value, ...)
%
% type is the observable:
% 'Nearest neighbor distance'
% 'Polarization'
% 'Speed'
% 'Distance covered'
% 'Coord pos'
% 'Absolute turn rate'
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
% 'Tail beat freq'
% 'Tail beat amplitude'
% 'Distance between ids'
% 'Transfer entropy'
%
% Other parameters that are needed based on the observable
%
% 'Fps' = frame rate
% 'FocalId' = id of the fish if you want to analyse a single fish. if not
% specified, group average is reported, specify two ids as [1 2], for xcorr 
% or transfer entropy, specify select ids for polarization, ANND, etc.
% 'NumTargets' = number of targets
% 'TrialTimeMin' = trial time in minutes. if this is a 1x2 vector then only
% the time between those two values in minutes is extracted
% 'RoiCm' = roi is a 1x2 vector e.g. [56 30]. a 1x4 vector assumes exact
% roi, e.g. [-27 -15 27 15] (xmin ymin xmax ymax]
% 'StimPosCm'= position of stimulus if only one value is specified then
% distance along the same is computed e.g. [10] or [5 10] 
% 'RegionBounds' = Boundaries of a region to compute number of entries;
% specify min and maximum e.g. [15, 28], or [-28, -15]
% 'TrajSec' = segments of trajectory in seconds that must be used to x corr
% 'TankPartitions' = number of partitions you want to divide the tank into
% 'Coord' = coordinate dimension along which the quantity must be computed
%
% Examples:
% 1) find the percentage time spent freezing by fish identified by id 1
% get_obs('Freezing', 'fps', 24, 'numTargets', 5,'TrialTimeMin', 5, 'Focalid', 1)
%
% 2) find the percentage time spent freezing by all fish
% get_obs('Freezing', 'fps', 24, 'numTargets', 5, 'RoiCm', [50 30], 'TrialTimeMin', 5)
% or a customized version with radius and deltawall
% get_obs('Freezing', 'fps', 60, 'numTargets', 1', 'TrialTimeMin', 1, 'RoiCm', [25.4 25.4], 'FreezingRadius', 2, 'deltaWall', 2)
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
%
% 6) find leadership time lag using xcorr
% get_obs('Leadership time lag xcorr', 'fps', 15, 'numTargets', 2,'TrialTimeMin', 5, 'FocalId', [1 2], 'TrajSec', 2);
% if 'FocalId' has a single id then the rest of the shoal is considered as
% the second id
%
% 7) find distance between two fish along the x-axis (Coord=1), if Coord
% not specified then distance is computed in the plane
% get_obs('Distance between ids', 'fps', 25, 'numTargets', 2, 'TrialTimeMin', 5, 'FocalId', [1 2], 'Coord', [1])
%
% 8) find the polarization between the targets
% get_obs('Polarization', 'fps', 25, 'numTargets', 2, 'TrialTimeMin', 5)
%
% 9) find the polarization between select target ids [1 2 4]
% get_obs('Polarization', 'fps', 25, 'FocalIds', [1 2 4], 'TrialTimeMin', 5)
%
% 10) find the information flow b/w two targets using transfer entropy
% get_obs('Transfer entropy', 'fps', 25, 'numTargets', 2, 'TrialTimeMin', 5, 'FocalId', [1 2], 'Coord', 2)
% In the above e.g., the information flow is measured from target 1 to 2;
% to measure the information flow from 2 to 1, FocalId should be [2 1]
% You can also compute transfer entropy with different values of
% downsampling interval (default 30) and number of bins (default=7) as
% follows
% get_obs('Transfer entropy', 'fps', 25, 'numTargets', 2, 'TrialTimeMin', 5, 'FocalId', [1 2], 'TEDS', 50, 'TENB', 7, 'Coord', 2)
% 
% 11) find tail beat frequency and tail beat amplitude. they both depend on
% the fish size in pixels and the length of the trajectory that you use to run an fft
% on. Hence, the command is
% get_obs('Tail beat freq', 'fps', 60, 'numTargets', 1, 'TrialtimeMin', 2, 'FishLengthPix', 40, 'TrajSec', 2*60)
% get_obs('Tail beat amp', 'fps', 60, 'numTargets', 1, 'TrialtimeMin', 2, 'FishLengthPix', 40, 'TrajSec', 2*60)
%


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
p.addParamValue('TrajSec', 0);
p.addParamValue('TankPartitions', 3);
p.addParamValue('FreezingRadius', 2);
p.addParamValue('Coord', 0);
p.addParamValue('TEDS', 30);
p.addParamValue('TENB', 7);
p.addParamValue('deltaWall', 3);
p.addParamValue('FishLengthPix', 0);


p.parse(varargin{:});

fid=p.Results.FocalId;
fps=p.Results.Fps;
nt=p.Results.NumTargets;
nmin=p.Results.TrialTimeMin;
roi=p.Results.RoiCm;
stimpos=p.Results.StimPosCm;
regionbounds=p.Results.RegionBounds;
trajsec=p.Results.TrajSec;
tparts=p.Results.TankPartitions;
frad=p.Results.FreezingRadius;
coord=p.Results.Coord;
ds=p.Results.TEDS;
nb=p.Results.TENB;
dwall=p.Results.deltaWall;
fpix=p.Results.FishLengthPix;

%%%% processing
if numel(find(fid~=0))>1
    nt=numel(fid);
end

if fid(1) && numel(find(fid~=0))==1 && nt ~=1 && ~strcmp(type, 'Leadership time lag xcorr') && ...
                       ~strcmp(type, 'Distance between ids') && ...
                       ~strcmp(type, 'Transfer entropy')
    warning('MATLAB:numtargets','Number of expected targets should be 1 for a focal fish');
end


% validating
if ~fps
    error('Fps must be specified');
end

if ~nmin
    error('TrialTimeMin must be specified');
end

if ~fpix && (strcmp('Tail beat freq', type) || strcmp('Tail beat amp', type))
    error('Tail beat frequency and amplitude need fish length in pixels. See help.');
end

if ~sum(roi) && (strcmp('Freezing', type) || strcmp('Thrashing', type) || strcmp('Swimming', type))
    error('RoiCm must be specified to properly classify locomotory behavior');
end

if isnan(stimpos) && strcmp('Distance from stimulus', type)
    error('StimPosCm must be specified for Distance from stimulus')
end

if numel(fid)~=2 && strcmp('Distance between ids', type)
    error('Exactly two ids must be specified for distance between ids')
end

if isnan(nt) && ~fid(1) && strcmp('Nearest Neighbor distance', type)
    error('NumTargets must be specified for Nearest Neighbor Distance if no ids are specified')
end

if numel(fid)~=2 && (strcmp('Transfer entropy', type) || strcmp('Leadership time lag xcorr', type))
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
        get_data=@(PathName, FileName) get_annd(PathName, FileName, nt, nfrm, fps, fid);
    case 'Polarization'
        bin_centers=0:0.05:1; units=''; acr='pol'; 
        get_data=@(PathName, FileName) get_pol(PathName, FileName, nt, nfrm, fps, fid);
    case 'Speed'
        bin_centers=0:1:10; units='(cm/s)'; acr='speed'; 
        get_data=@(PathName, FileName) get_speed(PathName, FileName, nt, nfrm, fps, fid);    
    case 'Distance covered'
        bin_centers=0:1:100; units='(m)'; acr='dtcm'; 
        get_data=@(PathName, FileName) get_path_travelled(PathName, FileName, nt, nfrm, fps, fid);   
    case 'Coord position'
        bin_centers=0:5:100; units='(cm)'; acr='pos'; 
        get_data=@(PathName, FileName) get_coord_pos(PathName, FileName, nt, nfrm, fps, coord, fid);    
    case 'Absolute turn rate'
        bin_centers=0:5:300; units='(degree/s)'; acr='turnrate'; 
        get_data=@(PathName, FileName) get_turn_rate(PathName, FileName, nt, nfrm, fid, fps);        
    case 'Freezing'
        bin_centers=0:.1:1; units='(%)'; acr='ptf'; acc_t=120; tsec=2;
        get_data=@(PathName, FileName) get_freeze(PathName, FileName, nt, nfrm, fid, fps,...
                                                    tsec, acc_t, roi, dwall, frad); 
    case 'Thrashing'
        bin_centers=0:.1:1; units='(%)'; acr='ptt'; acc_t=120; tsec=2;
        get_data=@(PathName, FileName) get_thrash(PathName, FileName, nt, nfrm, fid, fps,...
                                                    tsec, acc_t, roi, dwall, frad);
    case 'Swimming'
        bin_centers=0:.1:1; units='(%)'; acr='ptt'; acc_t=120; tsec=2;
        get_data=@(PathName, FileName) get_swimm(PathName, FileName, nt, nfrm, fid, fps,...
                                                    tsec, acc_t, roi, dwall, frad);                                             
    case 'Excursions'
        bin_centers=0:.1:1; units='(%)'; acr='excur'; 
        get_data=@(PathName, FileName) get_excursions(PathName, FileName, nt, nfrm, fps, fid);    
    case 'Leadership time lag xcorr'
        bin_centers=-.5:.1:.5; units='(sec)'; acr='timelag';
        get_data=@(Pathname, FileName) get_leadership_xcorr(PathName, FileName, nt, nfrm, fps, trajsec, fid);
    case 'Distance between ids'
        bin_centers=0:5:100; units='(cm)'; acr='distbwids';
        get_data=@(Pathname, FileName) get_distbwids(PathName, FileName, nt, nfrm, fps, fid, coord);      
    case 'Time spent left'
        bin_centers=0:.1:1; units='(%)'; acr='tsl'; 
        get_data=@(PathName, FileName) get_tsl(PathName, FileName, nt, nfrm, fps, fid); 
    case 'Time spent right'
        bin_centers=0:.1:1; units='(%)'; acr='tsr'; 
        get_data=@(PathName, FileName) get_tsr(PathName, FileName, nt, nfrm, fps, fid);    
    case 'Acceleration'
        bin_centers=0:1:20; units='(cm/s^2)'; acr='acc'; 
        get_data=@(PathName, FileName) get_accel(PathName, FileName,nt, nfrm, fid, fps);
    case 'Distance from stimulus'
        bin_centers=0:2:100; units='(cm)'; acr='d2s'; 
        get_data=@(PathName, FileName) get_dist2stim(PathName, FileName, nt, nfrm, fps, fid, stimpos);  
    case 'Time near stimulus'
        bin_centers=0:.1:1; units='(%)'; acr='tsns'; 
        get_data=@(PathName, FileName) get_time_near_stim(PathName, FileName, nt, nfrm, fps, fid, stimpos, roi, tparts);  
    case 'Time far from stimulus'
        bin_centers=0:.1:1; units='(%)'; acr='tsfs'; 
        get_data=@(PathName, FileName) get_time_far_stim(PathName, FileName, nt, nfrm, fps, fid, stimpos, roi, tparts);
    case 'Number of entries into a region'
        bin_centers=0:2:50; units=''; acr='freqr'; 
        get_data=@(PathName, FileName) get_freq_into_region(PathName, FileName, nt, nfrm, fps, fid, regionbounds);
    case 'Tail beat freq'
        bin_centers=0:1:10; units='(Hz)'; acr='tbf'; 
        get_data=@(PathName, FileName) get_tbf(PathName, FileName, nt, nfrm, fps, fid, trajsec, fpix);
    case 'Tail beat amp'
        bin_centers=0:1:1; units='(pixels)'; acr='tba'; 
        get_data=@(PathName, FileName) get_tba(PathName, FileName, nt, nfrm, fps, fid, trajsec, fpix);    
    case 'Transfer entropy'
        bin_centers=0:.2:2; units='(bits)'; acr='te12'; 
        get_data=@(PathName, FileName) get_trent(PathName, FileName, nt, nfrm, fps, coord, ds, nb, fid);
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% anal='Variation in '; 
% anal='Average ';

% figure(1); gcf; clf;
% ha=tight_subplot(1, numel(cselect), [0.025, 0.01], [0.2, 0.025], [0.05 .01]);


% ad_test_1way=nan(10, max(cselect));
% ad_test_2way=nan(10, 5, max(cselect));

[data permin]=get_data(PathName, FileName);

for jj=1:size(FileName,2)
    tmp=textscan(FileName{jj}, '%s', 'Delimiter', '_');
    if size(tmp{1},1) > 1
        cond{jj}=tmp{1}{2};
    else
        cond{jj}=tmp{1}{1};
    end
end
fprintf('\nAverage *%s* %s [mean value, 1st min, 2nd min, 3rd min, ... ]\n', type, units);
stats=nanmean(data,2);
% fprintf('FileName\t\t Avg. Value \t\t P1 \t P2 \t P3 \t P4 \t P5\n', type, units);
fprintf('------------------------------------------------\n');
for jj=1:size(FileName,2)
    fprintf('%s\t %.3f\t', cond{jj}, stats(jj));
    for ii=1:size(permin,2)
        fprintf('\t %.3f', permin(jj,ii));
    end
    fprintf('\n');
%     fprintf('%s\t\t %.3f \t\t %.3f \t %.3f \t %.3f \t %.3f \t %.3f\n', FileName{1,jj}, stats(jj), permin(jj,1), permin(jj,2), ...
%                                                                         permin(jj,3), permin(jj,4), permin(jj,5));
end

uav=~isnan(data);
fprintf('%% data used = %.2f +/- %.2f\n', mean(mean(uav,2))*100, std(mean(uav,2)*100));


plot(data');
xlabel('time (min)');
ylabel([type, units]);
[savedata, path, resp]=uiputfile(['dat0_', type, '.csv'], 'Save data file');
if resp, csvwrite([path, '/', savedata], data); end
if resp, csvwrite([path, '/', 'dat1_', type, '.csv'], [grp2idx(cond), stats, permin]); end


%% stats
fprintf('Statistical analysis: assumes that the data files are named as "X00_<condition name>_*"\n');
resp=input('Do you want to run One-way ANOVA stats on this data? []=yes, other=no: ');

if isempty(resp)

    fprintf('The stats are run with the following conditions as per the naming of files: \n')
    cc=unique(cond);
    condnum=grp2idx(cond);
    ccnum=unique(condnum);
    for ii=1:size(cc,2)
        fprintf('%s (%d files)\n',cc{ii}, sum(condnum==ii));
    end
    figure(2); gcf; clf;
    subplot(1,2,1); gca;
    set(gca, 'fontsize', 16);
    mystat('anovan', stats, cond, [type, units], 2);
    
    subplot(1,2,2); gca;
    for ii=ccnum'
        cdata=permin(ccnum==ii,:);
        mu=nanmean(cdata,1);
        st=nanstd(cdata,[],1);
        errorbar(mu, st, '-o', 'color', rand(1,3)); hold on
    end
    set(gca, 'fontsize', 16);
    legend(cc);
    xlabel('time');
    ylabel([type, units]);
end

function flag=filter_nt(rf, nt)
flag=double(size(rf,2)==nt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nearest Neighbor Distance
function [dat permin]=get_annd(PathName, FileName, nt, nfrm, fps, id)

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

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
if ~mod(nmin,1)
    permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
    if size(permin,2)==1, permin=permin'; end
else
    permin=[];
end
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


%% Speed
function [dat permin]=get_speed(PathName, FileName, nt, nfrm, fps, id)

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

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end

%% Distance covered
function [dat permin]=get_path_travelled(PathName, FileName, nt, nfrm, fps, id)

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    % to remove noise this is smoothened more
    for jj=1:size(data.X,1)
        data.X(jj,:)=sma(data.X(jj,:), ceil(fps));
    end

    if id
        r=getind(Xi.nX, 1, id, 1:Xi.nX,1);
        data.X=data.X(r,:);
        nt=1;
    end
     
     
%     for jj=1:size(data.X,1)
%         data.X(jj,:)=sma(data.X(jj,:), 15);
%     end
    
    
    r1=data.X(1:Xi.nX:end,:);
    r2=data.X(2:Xi.nX:end,:);
    dr1=diff(r1,1,2);
    dr2=diff(r2,1,2);
   
    nz=max(sum(r1(:,1:size(dr1,2))~=0,1),1);
   
    %dat(ii,1:size(dr1,2))=cumsum(sum(sqrt(dr1.^2+dr2.^2),1),2)./nz;
    dat(ii,1:size(dr1,2))=sum(sum(sqrt(dr1.^2+dr2.^2),1)./nz);

end

dat=dat/100;

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end


%% Coord position
function [dat permin]=get_coord_pos(PathName, FileName, nt, nfrm, fps, coord, id)

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    % to remove noise this is smoothened more
    for jj=1:size(data.X,1)
        data.X(jj,:)=sma(data.X(jj,:), ceil(fps));
    end

     if id
        r=getind(Xi.nX, 1, id, 1:Xi.nX,1);
        data.X=data.X(r,:);
        nt=1;
     end
    
    
    r1=data.X(coord:Xi.nX:end,:);
    %dat(ii,1:size(dr1,2))=cumsum(sum(sqrt(dr1.^2+dr2.^2),1),2)./nz;
    dat(ii,1:size(r1,2))=r1;

end

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end


% Swimming straight
% function [dat permin]=get_swim_straight(PathName, FileName, nt, nfrm, id, fps)
% 
% if numel(nfrm)==1
%     dat=zeros(size(FileName,2), nfrm);
% else
%     dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
% end
% for ii =1:size(FileName,2)
%     [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
%  
%     data.X= clean_and_smooth(data.X,Xi,nt);
% 
%      if id
%         r=getind(Xi.nX, 1, id, 1:Xi.nX,1);
%         data.X=data.X(r,:);
%         nt=1;
%      end
%     
%     
%     v1=data.X(3:Xi.nX:end,:);
%     v2=data.X(4:Xi.nX:end,:);
%    
%     % get current timestep and previous timestep
%     v1k=v1(:,1:end-1);
%     v1km1=v1(:,2:end);
%     
%     v2k=v2(:,1:end-1);
%     v2km1=v2(:,2:end);
%     
%     % normalize
%     v1kn=v1k./sqrt(v1k.^2+v2k.^2);
%     v2kn=v2k./sqrt(v1k.^2+v2k.^2);
%     
%     v1km1n=v1km1./sqrt(v1km1.^2+v2km1.^2);
%     v2km1n=v2km1./sqrt(v1km1.^2+v2km1.^2);
%     
%     dp=v1kn.*v1km1n+v2kn.*v2km1n;
%     
%     nz=sum(v1k~=0,1);
%    
%     tr=acos(dp)*fps*180/pi;
%     tr(v1k==0 & v2k==0)=0;
%     dat(ii,1:size(tr,2))=sum(real(tr),1)./nz;
%     dat(ii,nz~=nt)=nan;
% end
% 
% if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
% nmin=nfrm/fps/60;
% permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
% if size(permin,2)==1, permin=permin'; end


%% Absolute turn rate
function [dat permin]=get_turn_rate(PathName, FileName, nt, nfrm, id, fps)

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    % to remove noise this is smoothened more
%     for jj=1:size(data.X,1)
%         data.X(jj,:)=sma(data.X(jj,:), ceil(fps)/4);
%     end

    if id
        r=getind(Xi.nX, 1, id, 1:Xi.nX,1);
        data.X=data.X(r,:);
        nt=1;
    end
    
    v1=data.X(3:Xi.nX:end,:);
    v2=data.X(4:Xi.nX:end,:);
   

    % add heading info (if available)
    if Xi.nX ==10
        v1=data.X(9:Xi.nX:end,:);
        v2=data.X(10:Xi.nX:end,:);
    end
    
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
   
    tr=(acosd(dp)*fps);
    tr(v1k==0 & v2k==0)=0;
    dat(ii,1:size(tr,2))=sum(real(tr),1)./nz;
    dat(ii,nz~=nt)=nan;
end

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end

%% Polarization
function [dat permin]=get_pol(PathName, FileName, nt, nfrm, fps, id)

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

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;

if ~mod(nmin,1)
    permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
    if size(permin,2)==1, permin=permin'; end
else
    permin=[];
end


function P=pol(v)
% function pol(v)
%
% v is a d x n matrix with d dimensions and n targets
[d n]=size(v);

vhat = v./(ones(d,1)*sqrt(sum(v.^2)));

P=1/n*norm(sum(vhat,2));



%% Freezing 
function [dat permin]=get_freeze(PathName, FileName, nt, nfrm, id, fps, tsec, acc_t, roi, dwall, rad)

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
    fr1=locomotory_behavior(data.X, Xi.nX, fps, tsec, acc_t, roi, dwall, rad, Xi); 
    
    dat(ii,1:size(fr1,2))=fr1;
end

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end

% Thrashing 
function [dat permin]=get_thrash(PathName, FileName, nt, nfrm, id, fps, tsec, acc_t, roi, dwall, rad)

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
    [fr1 tr1]=locomotory_behavior(data.X, Xi.nX, fps, tsec, acc_t, roi, dwall, rad, Xi); 

    dat(ii,1:size(tr1,2))=tr1;
end

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end

% Swimming 
function [dat permin]=get_swimm(PathName, FileName, nt, nfrm, id, fps, tsec, acc_t, roi, dwall, rad)

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
    [fr1 tr1 sw1]=locomotory_behavior(data.X, Xi.nX, fps, tsec, acc_t, roi, dwall, rad, Xi); 

    dat(ii,1:size(sw1,2))=sw1;
end

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end


function [freeze thrash swimm]=locomotory_behavior(X, nX, fps, tsec, acc_t, roi, dwall, dist_t, Xi)
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

% add heading info
if Xi.nX ==6
    h1=v1;
    h2=v2;
else
    h1=X(9:Xi.nX:end,:);
    h2=X(10:Xi.nX:end,:);
end
heading=angle(h1+1i*h2);

% for later 
turnRate=diff(heading, 1, 2)*fps;

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
        if max(dist(1,:))<dist_t && ~sum(dist==0) % if the fish stayed within a ball of radius dist_t
            if numel(roi)==2
                if sum(accel(ff,ii:ii+n-1)) > acc_t && ... % if the fish is accelerating above a value 
                    (sum(abs(rt(1,:))>roi(1)/2-dwall) || ... % if the fish is near the wall, i.e., 
                    sum(abs(rt(2,:))>roi(2)/2-dwall))        % the fish position is within a delta of the wall
                    thrash(ii:ii+n-1)=1;
                else
                    freeze(ii:ii+n-1)=1;
                end
            elseif numel(roi)==4
                if sum(accel(ff,ii:ii+n-1)) > acc_t && ...
                    (sum(rt(1,:)<roi(1)+dwall) || ...
                     sum(rt(2,:)<roi(2)+dwall) || ...
                     sum(rt(1,:)>roi(3)-dwall) || ...
                     sum(rt(2,:)>roi(4)-dwall))
                    thrash(ii:ii+n-1)=1;
                else
                    freeze(ii:ii+n-1)=1;
                end
            % this new is based on Violet and Tiziana's analysis that fish can also be considered
            % thrashing if it is at an angle of 20 degrees or more to the
            % wall            
%             elseif (sum(abs(rt(1,:))>roi(1)/2-dwall) && acosd(abs(h1))>20 && acosd(abs(h1))<160) || ...
%                    (sum(abs(rt(2,:))>roi(2)/2-dwall) && acosd(abs(h1))<70)
%                 thrash(ii:ii+n-1)=1;
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


%% Excursions
function [dat permin]=get_excursions(PathName, FileName, nt, nfrm, fps, id)

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end

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

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end

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



%% leadership using xcorr
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
        if numel(xcorrids)==2
            r2e=getind(Xi.nX, k, xcorrids(2), 3:4,1);
            v2e=normv(data.X(r2e,k:k+trajsec));
            v1us=v2e(1,:);
            v2us=v2e(2,:);
        else
            % [TODO] fix this so that any xcorrid(1) may be picked!
            v1=data.X(3:Xi.nX:end,k:k+trajsec);
            v2=data.X(4:Xi.nX:end,k:k+trajsec);
            nz=sum(v1~=0,1);

            v1us=sum(v1(2:end,:),1)./(nz-1);
            v2us=sum(v2(2:end,:),1)./(nz-1);
            
        end

        
        
        dir_us=atan2(v2us, v1us);
        
%         [xx lags]=xcorr(dir_e, dir_us, 'coeff');
        [xx lags]=xcorr(dir_e, dir_us);
%         [xx lags]=xcorr(vve(1,:), v1us);
        if ~isnan(max(xx))
            [val idx]=max(xx);
            tau=lags(idx)/fps;
        else
            tau=0;
        end
        
        
%         tau=lags(idx)/fps*val;
        dat(ii,k:k+trajsec)=tau;
    end
end

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end

%% Time spent left
function [dat permin]=get_tsl(PathName, FileName, nt, nfrm, fps, id)

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

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end

%% Time spent right
function [dat permin]=get_tsr(PathName, FileName, nt, nfrm, fps, id)

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

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end

%% Acceleration
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

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end

%% Distance from stimulus
function [dat permin]=get_dist2stim(PathName, FileName, nt, nfrm, fps, id, stimpos)

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

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end



%% Distance between ids
function [dat permin]=get_distbwids(PathName, FileName, nt, nfrm, fps, ids, coord)

% if coord is not specified then compute on both dimensions
if ~coord, coord=1:2; end;

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end

for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);

    data.X= clean_and_smooth(data.X,Xi,nt);

    
    r1=getind(Xi.nX, 1, ids(1), 1:2,1);
    nzidx1=find(data.X(r1(1),:)~=0);
    r2=getind(Xi.nX, 1, ids(2), 1:2,1);
    nzidx2=find(data.X(r2(1),:)~=0);

    commonidx=intersect(nzidx1,nzidx2);

    dat(ii,commonidx)=sqrt(sum((data.X(r1(coord),commonidx)-data.X(r2(coord),commonidx)).^2,1));
    if numel(commonidx)<0.9*size(r1,2)
        fprintf('[?] the two fish are not tracked for the same frames...\n');
    end
end

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end

%% Time near stimulus
function [dat permin]=get_time_near_stim(PathName, FileName, nt, nfrm, fps,id, stimpos, roi, tparts)

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

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end

%% Time far from stimulus
function [dat permin]=get_time_far_stim(PathName, FileName, nt, nfrm, fps, id, stimpos, roi, tparts)

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

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end

%% Entries into a region
function [dat permin]=get_freq_into_region(PathName, FileName, nt, nfrm, fps, id, regionbounds)

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

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(sum(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end

dat=sum(dat,2);

%% tail beat amplitude
function [dat permin]=get_tba(PathName, FileName, nt, nfrm, fps, id, trajsec, fpix)

debug=0;

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
    
   

    tba_1=nan(nt, size(data.X,2));
    for ff=1:nt
        [tbf, tba_1(ff,:)]=tbf_perfish(data.X, Xi, ff,fps, trajsec, fpix);
    end
    
    
%     for k=1:size(data.X,2)
%         [rf vf hf]=processk(data.X, Xi, k,0);
%         P(ii,k)=pol(hf(1:2,:));
%         if size(rf,2)~=nf, P(ii,k)=nan; end
%         if isnan(P(ii,k)), keyboard; end
%     end
    dat(ii,:)=nanmean(tba_1,1);
end
if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end


%% tail beat frequency
function [dat permin]=get_tbf(PathName, FileName, nt, nfrm, fps, id, trajsec, fpix)

debug=0;

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
    
    if debug
        figure(1); gcf; clf;
        xx=-100:1:0;
        for k=5:size(data.X,2)-5
            for ff=1:nt
                subplot(1,3,ff); gca; cla;
                r=getind(Xi.nX, 1, ff, 1:Xi.nX,1);
                for cc=1:5
                    if data.X(r(Xi.sh(1)),k+cc) ~=1 
                        plot(xx, data.X(r(Xi.sh(1)),k+cc)*xx.^2 + data.X(r(Xi.sh(2)),k+cc)*xx, 'Color', zeros(1,3)+cc/10);
                    end
                    hold on;
                end
                set(gca, 'ylim', [-100 100]);
                set(gca, 'xlim', [-100 0]);
            end
            xlabel(k);
            drawnow;
        end
    end

    tbf_1=nan(nt, size(data.X,2));
    for ff=1:nt
        tbf_1(ff,:)=tbf_perfish(data.X, Xi, ff,fps, trajsec, fpix);
    end
    
    
%     for k=1:size(data.X,2)
%         [rf vf hf]=processk(data.X, Xi, k,0);
%         P(ii,k)=pol(hf(1:2,:));
%         if size(rf,2)~=nf, P(ii,k)=nan; end
%         if isnan(P(ii,k)), keyboard; end
%     end
    dat(ii,:)=nanmean(tbf_1,1);
end
if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end


function [tbf_1 tba_1]=tbf_perfish(X, Xi, ff, fps, trajsec, fpix)

dbg=0; % use this to plot the procedure for tbf

%%%%
% Take the shape parameters which is simply the coefficients of a quadratic
% curve that fit the fish shape. Find the tip of the tail on this curve and
% use the time-series to find the frequency. In particular, use fft on
% sequences that are longer than 30 frames and where (a) the fish is
% oriented against the flow and (b) it is not close to the wall of the
% tanks. After detrending the data to remove stationarity effects, use 
% the frequency at the maximum power as the tail-beat frequency of a single
% fish
%

tbf_1=nan(1,size(X,2));
tba_1=nan(1,size(X,2));

% this value is in pixels on the screen; may change depending on the experimental setup
xx=-fpix/2; 

Fs=fps; % sampling frequency

r=getind(Xi.nX, 1, ff, 1:Xi.nX,1);
a=X(r(Xi.sh(1)),:); 
b=X(r(Xi.sh(2)),:);

% the furthest point on the tip of the tail
ey=a*xx^2 + b*xx;


% pick sections of 2*fps
for dd=1:trajsec*fps:numel(ey)-trajsec*fps
    idx=dd:dd+2*fps-1;
    ey1=ey(idx);

    ey1=detrend(ey1,0);
    freq=Fs./length(ey1)*(1:numel(ey1));
    mag=abs(fft(ey1))./length(ey1);
    [val idx2]=max(mag);

%     L=numel(ey1);
%     NFFT = 2^nextpow2(L); % number of points on the fft
%     Y = fft(ey1,NFFT)/L;
%     Fs=fps; % sampling frequency
%     f = Fs/2*linspace(0,1,NFFT/2+1);
% 
%     pwr=2*abs(Y(1:NFFT/2+1));
    if dbg
        figure(10); gcf; clf;
        xb=0:-1:xx;
        a1=a(idx);
        b1=b(idx);
        subplot(1,3,1); gca;
        yb=(xb.^2)'*a1 + (xb.^2)'*b1;
        for jj=1:numel(a1)
            plot(xb, yb(:,jj), 'color', ones(1,3)*(.35+.65*jj/numel(a1)), 'linewidth', 3);
            hold on;
        end
        axis([xx 0 xx -xx]);
        set(gca, 'fontsize', 24);
        xlabel('pixels'); ylabel('pixels');
        
        
        subplot(1,3,2); gca;
        plot((1:numel(a1))/fps, ey1, 'b', 'linewidth', 3);
        box off;
        set(gca, 'fontsize', 24);
        ylabel('relative tail-tip position')
        xlabel('time (s)');
        
        
        subplot(1,3,3); gca;
        plot(freq, mag); hold on;
%         plot(f,pwr, 'linewidth', 3)
        set(gca, 'fontsize', 24);
%         title('Single-Sided Amplitude Spectrum of y(t)')
        xlabel('Frequency (Hz)')
        ylabel('Magnitude')
        plot(freq(idx2), mag(idx2), 'bx', 'linewidth', 4, 'markersize', 10);
        title(sprintf('%.1f Hz, %.1f amplitude', freq(idx2), mag(idx2)*2));
        pwr(1:2)=0;
        drawnow;
    end

    tbf_1(idx)=freq(idx2);
    tba_1(idx)=mag(idx2)*2; % this is in pixels??
%     [val idx2]=max(pwr);
%     tbf_1(idx)=f(idx2);
%     tba_1(idx)=val;
end


% pick the largest
% [v idx]=max(nv);
% tbf_1=[];
% for dd=1:numel(idx)
%     ey=eyy(idx(dd)).v;
%     ey=detrend(ey,0);
% 
% 
%     L=numel(ey);
%     NFFT = 2^nextpow2(L);
%     Y = fft(ey,NFFT)/L;
%     Fs=fps;
%     f = Fs/2*linspace(0,1,NFFT/2+1);
% %             subplot(1,2,2); gca;
%     pwr=2*abs(Y(1:NFFT/2+1));
% %             plot(f,pwr)
% %             title('Single-Sided Amplitude Spectrum of y(t)')
% %             xlabel('Frequency (Hz)')
% %             ylabel('|Y(f)|')
% %             pwr(1:2)=0;
%     [val idx2]=max(pwr);
%     tbf_1=[tbf_1 f(idx2)];
% end


%% Transfer entropy (information flow)
function [dat permin]=get_trent(PathName, FileName, nt, nfrm, fps, coord, ds, nb, teids)

if numel(teids)~=2
    error('both fish ids must be specified..');
end

if numel(nfrm)==1
    dat=zeros(size(FileName,2), nfrm);
else
    dat=zeros(size(FileName,2), nfrm(2)-nfrm(1));
end

labels={'fish (X)', 'stimulus (Y)'};
for ii =1:size(FileName,2)
    [data.X Xi]=reformat_csv([PathName, FileName{ii}], nfrm);
 
    data.X= clean_and_smooth(data.X,Xi,nt);
    
    rX=getind(Xi.nX, 1, teids(1), coord,1);
    X=data.X(rX,:);
    
    rY=getind(Xi.nX, 1, teids(2), coord,1);
    Y=data.X(rY,:);
    
    TE=transfer_entropy(X, Y, 1, 1, ds, nb, labels, 0); % show=0
    dat(ii,:)=TE.XY;
end

if numel(nfrm)==2, nfrm=nfrm(2)-nfrm(1); end
nmin=nfrm/fps/60;
permin=squeeze(nanmean(reshape(dat, [size(dat,1), size(dat,2)/nmin, nmin]), 2));
if size(permin,2)==1, permin=permin'; end

function TE=transfer_entropy(X,Y,k,l, ds, nb, labels, show)

% Schreiber, T., 2000. Measuring information transfer. Physical Review Letters, 85(2), pp.461?464.
% TE=p(X[n], X[n-1], ..., X[n-k], Y[n-1], Y[n-2],
% Y[n-l])log((p(X[n]|X[n-1], ..., X[n-k], Y[n-1], Y[n-2], ...,
% Y[n-l])/p(X[n]|X[n-1], ..., X[n-k]))
%
% k=1, l=1
%
% TE=p(X[n], X[n-1], Y[n-1])log((p(X[n]|X[n-1], Y[n-1])/p(X[n]|X[n-1]))
%
%

%%% toy problems
toy=0;
if nargin ==0
    X=1; toy=1;
end
if nargin==1
    toy=1;
end
if toy
    switch X
        case 1 % 
            t=0:.05:8*pi;
            X=sin(t)+randn(1,numel(t))*.05;
            Y=randn(1,numel(t))*.05;
            labels={'X', 'Y'};
        case 2
            ds=1; nb=1;
            B1=load('../data/B1.dat');
            X=B1(2350:ds:3550,1)';
            Y=B1(2350:ds:3550,2)';
            labels={'heart rate (X)', 'breath rate (Y)'};
        case 3
    end
    k=1;
    show=1;
end

% normalize
X=normalize_data(X);
Y=normalize_data(Y);

if show
    gca;
    plot(X, 'r');
	hold on;
    plot(Y, 'b');
    set(gca, 'fontsize', 18);
    xlabel('time');
end
r=nb;
TE.YX=trent_raw2(X,Y,k,ds, r);
TE.XY=trent_raw2(Y,X,k,ds, r);
TE.r=r;

% r=0.1
% TE.YX=trent_ci(X,Y,k,ds,r);
% TE.XY=trent_ci(Y,X,k,ds,r);
% TE.r=r;

if show
    gca;
    title(sprintf('TE(X->Y)=%.3f bits\t\t TE(Y->X)=%.3f bits', TE.XY, TE.YX));
    legend(labels);
end

function X=normalize_data(X)

% normalize
X=X-min(X(:)); X=X/max(X(:)); X=X*2-1;


function TE=trent_raw2(X,Y,k,ds, r)

% downsample
X=X(1:ds:end);
Y=Y(1:ds:end);

% binning
% Z=[X(:), Y(:)]; % use both datasets
% Z=X(:);
% bins=linspace(min(Z(:)), max(Z(:)), r);
bins=linspace(-1, 1, r); % because the incoming data is normalized anyway
% fprintf('binends=[%.1f, %.1f]\n', bins(1), bins(end));
binwidth=bins(2)-bins(1);
bb=bins;
bins=[-inf bins(1)-binwidth/2 bins+binwidth/2];
bins(end)=bb(end);
nbins=numel(bins);

% initialize the data
pXXY=zeros(nbins-1, nbins-1, nbins-1);
N=numel(X);

% create the markov chain
Xchain=zeros(k+1, N-k); 
for jj=k+1:-1:1
    Xchain(k+2-jj,:)=X(jj:jj+N-k-1);
end
Y=Y(1:N-k);


% compute 3D histogram
% Y=sin(2*pi*(1:numel(Y))/4); % this matches
% Y=sin(2*pi*(1:numel(Y))/4); % this is lower freq
for ii=1:nbins-1
    for jj=1:nbins-1
        for kk=1:nbins-1
            % the same k is used since we are using a logical operator
            pXXY(ii,jj,kk)=sum(Xchain(1,:)>=bins(ii) & Xchain(1,:)<=bins(ii+1) & ...
                               Xchain(2,:)>=bins(jj) & Xchain(2,:)<=bins(jj+1) & ...
                               Y>=bins(kk) & Y<=bins(kk+1));
        end
    end
end

pXXY=pXXY(2:end, 2:end, 2:end);
pXXY=pXXY/sum(pXXY(:));

% p(X[k-1], Y[k-1])
pXY=zeros(nbins-1, nbins-1);
for jj=1:nbins-1
    for kk=1:nbins-1
        pXY(jj,kk)=sum(Xchain(2,:)>=bins(jj) & Xchain(2,:)<=bins(jj+1) & ...
                               Y>=bins(kk) & Y<=bins(kk+1));
    end
end
pXY=pXY(2:end,2:end);
pXY=pXY/sum(pXY(:));


% p(X[k], X[k-1])
pXX=zeros(nbins-1,nbins-1);
for ii=1:nbins-1
    for jj=1:nbins-1
        pXX(ii,jj)=sum(Xchain(1,:)>=bins(ii) & Xchain(1,:)<=bins(ii+1) & ...
                               Xchain(2,:)>=bins(jj) & Xchain(2,:)<=bins(jj+1));
    end
end
pXX=pXX(2:end,2:end);
pXX=pXX/sum(pXX(:));


% p(X[k-1])
pXm1=zeros(nbins-1,1);
for jj=1:nbins-1
    pXm1(jj)=sum(Xchain(2,:)>=bins(jj) & Xchain(2,:)<=bins(jj+1));
end
pXm1=pXm1(2:end);
pXm1=pXm1/sum(pXm1(:));


pXXY(isnan(pXXY))=0;
pXX(isnan(pXX))=0;
pXY(isnan(pXY))=0;
pXm1(isnan(pXm1))=0;

% TE=p(X[n], X[n-1], Y[n-1])log((p(X[n]|X[n-1], Y[n-1])/p(X[n]|X[n-1]))
TE=0;
for ii=1:nbins-2
    for jj=1:nbins-2
        for kk=1:nbins-2
            if pXXY(ii,jj,kk) && pXX(ii,jj) && pXY(jj,kk) && pXm1(jj)
                TE=TE+pXXY(ii,jj,kk)*log2(pXXY(ii,jj,kk)/((pXX(ii,jj)/pXm1(jj))*pXY(jj,kk)));
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%________________________________________________________________

function [Xh Xi]=reformat_csv(csvfile, nfrm)

Xhcsv=csvread(csvfile);

% clean the file to ensure that some noise is removed here only
max_t=numel(unique(Xhcsv(:,2)));
if max_t > 100
    fprintf('[W] found more than 100 targets! deleting short tracks... \n');
    ids=unique(Xhcsv(:,2));
    ids2=[];
    for jj=ids'
        if sum(Xhcsv(:,2)==jj,1) > nfrm/10
            ids2=[ids2, jj];
        end
    end
    Xhcsv=Xhcsv(ismember(Xhcsv(:,2), ids2),:);
    max_t=max(Xhcsv(:,2));
end

if size(Xhcsv,2)>8
    if sum(Xhcsv(:,10))
        Xi=strXi('shape');
    else
        Xi=strXi;
    end
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

win=3;
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
    rf=[X(Xi.nX*(id-1)+1,k)';X(Xi.nX*(id-1)+2,k)'];
    if Xi.nX==10
        vf=[X(Xi.nX*(id-1)+9,k)';X(Xi.nX*(id-1)+10,k)'];
    else
        vf=[X(Xi.nX*(id-1)+3,k)';X(Xi.nX*(id-1)+4,k)'];
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

