function [bg nZ X P datfile cm2pix, sz1, mesg, frame]=onestep(sides_cm, img_t,...
                                bin_t, area_t, fgislight, circ, smooth, bg, alpha, ...
                                blur, bb_blur, bb_size, X, P, fps, ...
                                datfile, getfrm, bgyes, trktype, ...
                                write_tracks, record_verify, split, shape, sz1, frame)                           

global roi_crop
global roi_cut
global view_type
global k
global live1
global rp_id
global zoom_val
global trklen

mesg=struct('txt', '');
 
sz=[mean(sz1(sz1~=0)), std(sz1(sz1~=0))];

calib=calib2d(roi_crop, sides_cm);
cm2pix=@(x,y) cm2pix1(x,y, calib);

if live1
    gca; cla;
end

img=getfrm(k);

nZ=0;

[fg bg]=preproc(img, fgislight, smooth, img_t, bin_t, bg, alpha, ...
                    roi_crop, circ, roi_cut, blur, bb_blur, area_t, split,...
                    shape, sz, bb_size, frame);
if live1
    if ~bgyes
        imshow(img);
    else
        imshow(uint8(fg*255));
    end
end

% if ~isempty(bb)
% rectangle('position', bb);
% end

%%%%% MARK

if view_type==2
    if live1
        hold on;
        if mod(k,10)
            if circ
                center=[roi_crop(2)+roi_crop(4)/2; roi_crop(1)+roi_crop(3)/2];
                rad=roi_crop(3:4)/2;
                th=-pi:.1:pi;
                plot(center(2)+rad(1)*cos(th), center(1)+rad(2)*sin(th), 'b--');
            else
                rectangle('position', roi_crop, 'linestyle', '--', 'edgecolor', 'b');
            end
        end
        
        if mod(k,20);
            nzi=find(roi_cut(:,1)~=0);
            for jj=nzi'
                rectangle('position', roi_cut(jj,:), 'linestyle', '--', ...
                            'edgecolor', 'w');
            end
        end
    end
    [Zk ns_area bb frame]=getZ(fg, area_t, split, shape, sz, frame);
    nZ=size(Zk,2);
    % ensure that we are getting good shape parameters
    if shape
        if ~isempty(Zk)
            if sum(Zk(4,:)~=0)<size(Zk,2)
                mesg.txt='Shape :(';
            end
        end
    end
    if live1
        if nZ
            plot(Zk(1,:), Zk(2,:), 'r*');
            pix=mean(Zk(1:2,:),2);
            reset_axes_lim(pix', img, zoom_val);
        end
    end
end

%%%%% TRACK
if view_type==3
    if live1, hold on; end
    
    %%% get measurements
    [Zk1 ns_area bb frame]=getZ(fg, area_t, split, shape, sz, frame);
    
    % ensure that we are getting good shape parameters
    if shape
        if ~isempty(Zk1)
            if sum(Zk1(4,:)~=0)<size(Zk1,2)
                mesg.txt='Shape :(';
            end
        end
    end
    if ~isempty(Zk1)
        Zk1=Zk1(:,Zk1(1,:)~=0);
    end

    %%% change trackers here...
    try
        if shape && trktype ~=2
            error('[!] if tracking shape tracker type should be 2');
        end
        
        if trktype==0
            [X P]=mttkf2d(X, P,  Zk1, 1/fps, calib, 0);
        elseif trktype==1
            [X P]=mttkf2d(X, P,  Zk1, 1/fps, calib, 1);
        elseif trktype==2
%             error('shape tracker under maintenance!');
%             [X P]=mttpf2d(X, P,  Zk1, 1/fps, calib);
            [X P]=mttkf2ds(X, P,  Zk1, 1/fps, calib, 0);
        elseif trktype==3
            [X P]=mtt2d(X, P,  Zk1, 1/fps, calib, 0);    
        end
        frame(k).X=X(X(:,1)==k, :);
        frame(k).P=P(P(:,1)==k, :);
        
    catch ME
        fprintf('[!] Error; saving workspace ....\n');
        wsp=sprintf('./err_%s.mat', datestr(now, 'yyyymmddTHHMMSS'));
        fprintf('%s\n', ME.message);
        fprintf('line %d\n', ME.stack(1).line);
        fprintf('file %s\n', ME.stack(1).file);
        fprintf('%d', X);
        save(wsp);
        throw(ME);
    end
    
    if live1
        quick_verify(k,X, cm2pix);       
        pos=mean(X(X(:,1)==k,3:4),1);
        [xp yp]=cm2pix(pos(:,1), pos(:,2));
        reset_axes_lim([xp yp], img, zoom_val);
    end
    
    if write_tracks
        nzi=find(datfile(:,1)~=0); % find non zero entries
        if isempty(nzi), nzi=1; end
        curr=X(:,1)==k; % find current time-step
        nzi=nzi(end);
        curr_overwrite=datfile(:,1)==k; % delete if present
        if sum(curr_overwrite)
            datfile(curr_overwrite,:)=[];
            nzi=find(datfile(:,1)~=0); 
            if ~isempty(nzi)
                nzi=nzi(end); % update nzi
            else
                nzi=0;
            end
        end
        datfile(nzi+1:nzi+sum(curr),1:size(X,2))=X(curr,:);
    end
end


%%%%% REPAIR
if view_type==4    
    if live1, hold on; end
    quick_verify(k,datfile, cm2pix);
    

    pos=datfile(datfile(:,1)==k & datfile(:,2)==rp_id,3:4);

    if isempty(pos) && rp_id
        k_rp=datfile(datfile(:,2)==rp_id & datfile(:,1)<k,1);
        if ~isempty(k_rp)
            pos=datfile(datfile(:,1)==max(k_rp) & datfile(:,2)==rp_id,3:4);
        else
            pos=mean(datfile(datfile(:,1)==k,3:4),1);
        end
    elseif ~rp_id
        pos=mean(datfile(datfile(:,1)==k,3:4),1);
    end

    [xp yp]=cm2pix(pos(:,1), pos(:,2));
    reset_axes_lim([xp yp], img, zoom_val);
end


%%%%% VERIFY
if view_type==5
    if live1, hold on; end
    if live1
        if record_verify
            img1=quick_verify(k,datfile, cm2pix, img);
            
            imwrite(img1, sprintf('/tmp/verify_%.6d.jpg', k), 'jpg');
        else
            quick_verify(k,datfile, cm2pix);
            pos=mean(datfile(datfile(:,1)==k,3:4),1);
            [xp yp]=cm2pix(pos(:,1), pos(:,2));
            reset_axes_lim([xp yp], img, zoom_val);
        end
    end   
end

%%%%%%%%%%%%%%
if nZ
if numel(ns_area) 
    idx=randperm(1000);
    if numel(ns_area) < 1000
        sz1(idx(1:numel(ns_area)))=ns_area;
    end
end
end

if live1
    drawnow;
end


function [Zk ns_area bb frame]=getZ(fg, area_t, split, shape, sz, frame)

global k

bin_labels=logical(fg);
% regionprops/contours
Zk=[]; bb=[]; ns_area=[];
stats=regionprops(bin_labels,'Area', 'Centroid', 'PixelList', 'BoundingBox');
if ~isempty(stats)
    stats=stats(cat(1,stats.Area) > area_t);
    
    ns_area=cat(1,stats.Area);
    
    %%%%%%%%%%%% don't run if there are too many targets... probably due to
    %%%%%%%%%%%% bad contrast and threshold
    if size(stats,1) < 100
        stats=occlusion_splitting(stats, split, sz);
    end
    
    if shape
        [stats pixw]=curve_fit(stats); 
        Zk=[cat(1,stats.Centroid)';
            cat(1,stats.Area)';
            cat(2, stats.p);
            cat(2, stats.hd);];
    else
        Zk=[cat(1,stats.Centroid)';
            cat(1,stats.Area)';];
    end
    if ~isempty(frame(k).Mk)
        Zk=frame(k).Mk;
    end
    frame(k).Zk=Zk;
    bb=cat(1,stats.BoundingBox);
end

function [fg bg]=preproc(img, fgislight, smooth, img_t, bin_t, ...
                        bg, bg_rupdate_alpha, roi_crop, circ, roi_cut, ...
                        blur, bb_blur, area_t, split, shape, sz, bb_size, frame)

[ih, iw, id]=size(img);
if id > 1, img=rgb2gray(img); end

if fgislight
    img=255-img;
end

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

bb=[];
if bb_blur(3)
    fg=double(fg);
    [Zk ns_area bb frame]=getZ(fg, area_t, split, shape, sz, frame);
    
    if ~isempty(Zk)
        [val idx]=find(Zk(3,:)>bb_size);
        for jj=idx
            bb_iter=ceil(bb(jj,:));
            ext=bb_blur(4);
            bb_iter=[max(bb_iter(1)-ext,1), max(bb_iter(2)-ext,1), ...
                     min(bb_iter(3)+2*ext, iw-bb_iter(1)-1), ...
                     min(bb_iter(4)+2*ext, ih-bb_iter(2)-1)];
    %         try
            bb_fg=fg(bb_iter(2):bb_iter(2)+bb_iter(4), bb_iter(1):bb_iter(1)+bb_iter(3));
    %         catch
    %             keyboard
    %         end
            bb_fg=filter2(fspecial('gaussian', bb_blur(1:2), ...
                            bb_blur(3)), bb_fg);
            fg(bb_iter(2):bb_iter(2)+bb_iter(4), bb_iter(1):bb_iter(1)+bb_iter(3))=bb_fg;     
        end
    end
end


function [u v]=cm2pix1(x,y, calib)

u= x/calib.pix2cm(1)+calib.center(1);
v= y/calib.pix2cm(2)+calib.center(2);


function reset_axes_lim(pix, img, zoom_val)
global rp_handle

if zoom_val
    rp_handle=plot(pix(1), pix(2), 'cs', 'linewidth', 1);
    box_w=size(img,2)/zoom_val;
    xmin=max(1,ceil(pix(1)-box_w/2));
    xmax=min(floor(pix(1)+box_w/2), size(img,2));
    xlim=[xmin, xmax];

    box_h=size(img,1)/zoom_val;
    ymin=max(1,ceil(pix(2)-box_h/2));
    ymax=min(floor(pix(2)+box_h/2), size(img,1));
    ylim=[ymin, ymax];
else
    rp_handle=plot(pix(1), pix(2), 'cs', 'linewidth', 1);  
    xlim=[1, size(img,2)];
    ylim=[1, size(img,1)];
end

set(gca, 'xlim',xlim);
set(gca, 'ylim',ylim);

%%%% splitting occlusions
function stats = occlusion_splitting(stats, split, sz)
global debug

if isnan(sz(1))
    aa=cat(1,stats.Area);
    sz=[mean(aa) std(aa)];
end

ss=1;
if split 
    nst=size(stats,1);
    didx=[];
    if debug, clr=rand(10,3); end
    for ii=1:nst
        if stats(ii).Area > sz(1)+2*sz(2) % more than 95% 
            nf=round(stats(ii).Area/sz(1)); % number of fish in this occlusion
            if nf > 1
                didx=[didx, ii];
                px=stats(ii).PixelList; % get the pixels
                labels=emgm(px', nf); % split the occlusion into the number of fish
                for jj=1:nf
                    px1=px(labels==jj, :);
                    if ~isempty(px1)
                        stats(nst+ss,1).Centroid=mean(px1,1);
                        stats(nst+ss,1).Area=size(px1,1);
                        stats(nst+ss,1).PixelList=px1;
                        if debug
                            plot(px1(:,1), px1(:,2), '.', 'color', clr(jj,:), 'markersize', 1);
                            hold on;
                        end
                        ss=ss+1;
                    end
                end
            end
        end
    end
    stats(didx)=[];
end


function [stats pixw_out]=curve_fit(stats)
global debug
% just add two fields if no blobs are seen
if isempty(stats)
    [stats(:).p]=deal();
    [stats(:).hd]=deal();
end
pixw_out(1:size(stats,1))=struct('pts', []);
for jj=1:size(stats,1)
    % select the body pixels
    pix=stats(jj).PixelList;
    if ~isempty(pix) && size(pix,1) > 20 % the number of pixels must be enough to fit the curve
        % update the centroid to fall on the actual fish body
        old_centr=stats(jj).Centroid';
        [val idx]=min(sum((pix'-old_centr*ones(1,size(pix,1))).^2));
        body_pix=pix(idx,:)';
        shift_vec=body_pix-old_centr;
        if norm(shift_vec) > 1 % shift a bit more if the distance was more than 1
            stats(jj).Centroid=(old_centr+shift_vec*2)';
        end
%         old_centr-stats(jj).Centroid'
        % choose no more than 100 
        np=min(100,size(pix,1));
        pix=pix(ceil(linspace(1, size(pix,1), np)), :);
        T=linspace(-pi,pi,180);
    %     [P T]=meshgrid(1:size(pix,1), linspace(1,2*pi, 360));
        % for all combinations of theta and head, get the best fit
        np=numel(T);
        err=9999999;
        for kk=1:np
    %         hh=pix(P(kk),:);
            hd=[cos(T(kk)); sin(T(kk))];
            wTb=[hd, [-hd(2), hd(1)]', stats(jj).Centroid'; 0 0 1];
            bTw=inv(wTb);
            pixb=tra2b(pix', bTw);
            xx=pixb(1,:)'; yy=pixb(2,:)';
            % y = ax^2+bx (fit a parabola to the fish shape)
            A=[xx.^2, xx];
            b=yy;
            p=(A'*A)\A'*b;
            err1=norm(A*p-b);
            if err1<err
                % if the number of pixels in +ve x is less than -ve x that
                % means that the heading is flipped. based on the
                % assumption that the number of pixels are always more
                % towards the head
                if sum(pixb(1,:)>0) < sum(pixb(1,:)<0)
                    hd=-hd;
                    p(1)=-p(1);
                end
                err=err1;
                stats(jj).p=p;
                stats(jj).hd=hd;
            end
        end
    else
        stats(jj).p=zeros(2,1);
        stats(jj).hd=[1 0]';
    end
    if debug
        if ~isempty(stats)
            wTb=[stats(jj).hd, [-stats(jj).hd(2), stats(jj).hd(1)]', stats(jj).Centroid'; 0 0 1];
            pixb=tra2b(pix', inv(wTb));
            xx=pixb(1,:)';
            A=[xx.^2, xx];
            yy=A*stats(jj).p;
            pixw=tra2b([xx,yy]', wTb);
            quiver(stats(jj).Centroid(1), stats(jj).Centroid(2), ...
                            stats(jj).hd(1), stats(jj).hd(2), 10, 'r');
            plot(pixw(1,:), pixw(2,:), 'g.', 'markersize', 4);
            pixw_out(jj).pts=pixw;
        end
    end
end

