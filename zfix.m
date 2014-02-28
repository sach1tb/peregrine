function zfix(frmloc)
%function zfix(frmloc)
% frmloc is the location of frames
%
% IMPORTANT: this script expects that preferences have been set for the
% dataset using the peregrine GUI. Therefore, please run and set proper
% thresholds before running this script.

frmloc=[frmloc, '/'];


id='*.jpg'; % if you pass a direct argument

frmlist = dir([frmloc, id]);
getimage=@(k) imread([frmloc, frmlist(k).name]);
nfrm=size(frmlist,1);

suff=frmlist(1).name(1:end-6);
if exist([frmloc, 'config_', suff, '.csv'], 'file')
    conf=csvread([frmloc, 'config_', suff, '.csv']);
else
    error('check config');
end

bglist=dir([frmloc, 'background_*.*']);
if ~isempty(bglist)
    bg=imread([frmloc, bglist(1).name]);
else
    img=imread([frmloc, frmlist(1).name]);
    bg=rgb2gray(img);
end


nt=conf(1,1);
smooth=conf(14,:);
blur=conf(15,:);
bb_blur=conf(16,:);

roi_crop=conf(3,:);
            roi_cut=conf(4:13,:); 
roi_cut=roi_cut(roi_cut(:,1)~=0,:);

img_t=conf(2,1);
area_t=conf(2,2);
bin_t=conf(2,3);
circ=conf(18,4);
bg_rupdate_alpha=conf(19,2);


zfile=[frmloc, 'Z01_', suff, '.csv'];
if ~exist(zfile, 'file')
    Z1=zeros(3*nfrm,nt);
    csvwrite(zfile, Z1);
else
    Z1=csvread(zfile);
end

warning off;

for k=1:nfrm
    Zk1=Z1(3*k-2:3*k,:);
    if sum(Zk1(3,:)>0) < nt
        %if ~sum(Zk1(3,:)>200)

        img=imread([frmloc, frmlist(k).name]);
        [ih, iw, id]=size(img);
        if id > 1, img=rgb2gray(img); end


        if smooth(3), 
            filter2(fspecial('gaussian', smooth(1:2), smooth(3)), img); 
        end

        fg=img < img_t;
        if bin_t
            fg=fg & bg-img>bin_t;
            if bg_rupdate_alpha
                bg=bg*(1-bg_rupdate_alpha) + bg_rupdate_alpha*img;
            else
                bg=max(bg, img);
            end
        end
        if roi_crop(1)
            fg=setRoi(fg, roi_crop, roi_cut(roi_cut(:,1)~=0,:), circ);
        end


        if blur(3)
            fg=filter2(fspecial('gaussian', blur(1:2), ...
                            blur(3)), fg);
        %         fg=filter2(fspecial('disk', blur(3)), fg);
        end

        figure(1); gcf; clf;
        [h w]=size(img);
        imshow(img);
        
        hold on;
        set(gca, 'fontsize', 18);
        text(w*.75, h/10, ...
            sprintf('To add a point press [a] then mark\nTo delete a point press [d] then mark\nTo go to next frame press [space]\n'), ...
            'fontsize', 16);
        xlabel(k);
        

        Zk=getZ(fg, area_t);
        nZk=size(Zk,2);
        
        % delete if you see more than twice the number of points
        if nZk >=2*nt
            Zk=[]; nZk=0;
        end
        if nZk
            plot(Zk(1,:), Zk(2,:), 'r*');
        end

        ch=0;
        while ch~=32 % space bar to proceed but only if the number of targets are same as nt
            ch=getkey(1);
            if ch==97 % add a point
                [x y]=ginput(1);
                Zk(:,end+1)=[x;y;1];
                plot(x,y, 'r*');
            elseif ch==100 % delete a point
                [x y]=ginput(1);
                if nZk
                   dist=sum(([x;y]*ones(1, size(Zk,2)) - Zk(1:2,:)).^2);
                   [val idx]=min(dist);
                   plot(x,y, 'bo');  
                   Zk(:,idx)=[];
                end
            end
            Zk1=Zk;
        end

        if size(Zk1,2)==nt
            Z1(3*k-2:3*k,:)=Zk1;
            csvwrite([frmloc, 'Z01_', suff, '.csv'], Z1);
        else
            tmp=zeros(3,nt);
            tmp(:,1:size(Zk1,2))=Zk1;
            Z1(3*k-2:3*k,:)=tmp;
            csvwrite([frmloc, 'Z01_', suff, '.csv'], Z1);
        end
    end
    
end

warning on;


function [Zk bb]=getZ(fg, area_t)

bin_labels=logical(fg);
% regionprops/contours
Zk=[]; bb=[];
stats=regionprops(bin_labels,'Area', 'Centroid', 'BoundingBox');
if ~isempty(stats)
    stats=stats(cat(1,stats.Area) > area_t);
    Zk=[cat(1,stats.Centroid)';
        cat(1,stats.Area)';];
    bb=cat(1,stats.BoundingBox);
end


function fg = setRoi(fg, roi_crop, roi_cut, circ)

if circ % actually it's an ellipse
    center=[roi_crop(1)+roi_crop(3)/2; roi_crop(2)+roi_crop(4)/2];
    rad=roi_crop(3:4)/2;
    [H W]=meshgrid(1:size(fg,2), 1:size(fg,1));
    H=H-center(1);
    W=W-center(2);
    rad_p=W.^2/rad(2)^2 + H.^2/rad(1)^2;
    idx=rad_p>1;

    fg(idx)=0; %??
else
    % remove everything except roi_crop
    fg(1:roi_crop(2),:)=0;
    fg(roi_crop(2)+roi_crop(4):end,:)=0;
    fg(:,1:roi_crop(1))=0;
    fg(:,roi_crop(1)+roi_crop(3):end)=0;

end

for jj=1:size(roi_cut,1)
fg(roi_cut(jj,2):roi_cut(jj,2)+roi_cut(jj,4), ...
    roi_cut(jj,1):roi_cut(jj,1)+roi_cut(jj,3))=0;
end


function [ch, tim] = getkey(N,nonascii)

% GETKEY - get a keypress
%   CH = GETKEY waits for a single keypress and returns the ASCII code. It
%   accepts all ascii characters, including backspace (8), space (32),
%   enter (13), etc, that can be typed on the keyboard.
%   Non-ascii keys (ctrl, alt, ..) return a NaN. CH is a double.
%
%   CH = GETKEY(N) waits for N keypresses and returns their ASCII codes.
%   GETKEY(1) is the same as GETKEY without arguments.
%
%   GETKEY('non-ascii') or GETKEY(N,'non-ascii') uses non-documented
%   matlab features to return a string describing the key pressed.
%   In this way keys like ctrl, alt, tab etc. can also distinguished.
%   The return is a string (when N = 1) or a cell array of strings.
%
%   [CH,T] = GETKEY(...) also returns the time between the start of the
%   function and each keypress. This is, however, not that accurate.
%
%   This function is kind of a workaround for getch in C. It uses a modal,
%   but non-visible window, which does show up in the taskbar.
%   C-language keywords: KBHIT, KEYPRESS, GETKEY, GETCH
%
%   Examples:
%
%    fprintf('\nPress any key: ') ;
%    ch = getkey ;
%    fprintf('%c\n',ch) ;
%
%    fprintf('\nPress the Ctrl-key within 3 presses: ') ;
%    ch = getkey(3,'non-ascii')
%    if ismemmber('control', ch),
%      fprintf('OK\n') ;
%    else
%      fprintf(' ... wrong keys ...\n') ;
%    end
%
%  See also INPUT, UIWAIT
%           GETKEYWAIT (File Exchange)

% for Matlab 6.5 and higher
% version 2.0 (jun 2012)
% author : Jos van der Geest
% email  : jos@jasen.nl
%
% History
% 1.0 2005 - creation
% 1.1 dec 2006 - modified lay-out and help
% 1.2 apr 2009 - tested for more recent MatLab releases
% 1.3 jan 2012 - modified a few properties, included check is figure still
%            exists (after comment on FEX by Andrew).
% 2.0 jun 2012 - added functionality to accept multiple key presses

t00 = tic ; % start time of this function

% check the input arguments
error(nargchk(0,2,nargin))
switch nargin
    case 0
        nonascii = '' ;
        N = 1 ;
    case 1
        if ischar(N),
            nonascii = N ;
            N = 1 ;
        else
            nonascii = '' ;
        end
end

if numel(N) ~= 1 || ~isnumeric(N) || N < 1 || fix(N) ~= N
    error('N should be a positive integer scalar.') ;
end

% Determine the callback string to use
if strcmpi(nonascii,'non-ascii'),
    % non-ascii characters are accepted
    nonascii = true ;
    callstr = 'set(gcbf,''Userdata'',get(gcbf,''Currentkey'')) ; uiresume ' ;
elseif isempty(nonascii)
    nonascii = false ;
    % only standard ascii characters are accepted
    callstr = 'set(gcbf,''Userdata'',double(get(gcbf,''Currentcharacter''))) ; uiresume ' ;
else
    error('String argument should be the string ''non-ascii''') ;
end

% Set up the figure
% May be the position property  should be individually tweaked to avoid visibility
fh = figure(...
    'name','Press a key', ...
    'keypressfcn',callstr, ...
    'windowstyle','modal',...
    'numbertitle','off', ...
    'position',[0 0  1 1],...
    'userdata','timeout') ;
try
    ch = cell(1,N) ;
    tim = zeros(1,N) ;
    
    % loop to get N keypresses
    for k=1:N
        % Wait for something to happen, usually a key press so uiresume is
        % executed
        uiwait ;
        tim(k) = toc(t00) ; % get the time of the key press
        ch{k} = get(fh,'Userdata') ;  % and the key itself
        if isempty(ch{k}),
            if nonascii
                ch{k} = NaN ;
            else
                ch{k} = '' ;
            end
        end
    end
    if ~nonascii
        ch = [ch{:}] ;
    else
        if N==1
            ch = ch{1} ; % return as a string
        end
        % return as a cell array of strings
    end
 catch
    % Something went wrong, return empty matrices.
    ch = [] ;
    tim = [] ;
end

% clean up the figure, if it still exists
if ishandle(fh)
    delete(fh) ;
end
    