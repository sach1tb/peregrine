function make_movie(floc, imgid, datafile, configfile)
% function make_movie(floc, imgid, datafile, configfile)
%
% floc is the location where images are kept
% imgid is the image identifier e.g. CTRL_
% datafile is the X00_ csv datafile
% configfile is the config_ csv file
%
%

Xh=csvread(datafile);
flist=dir([floc, '*', imgid, '*']);

Xh=Xh(Xh(:,2)==1,:);
[val kidx]=sort(Xh(:,1));
Xh=Xh(kidx,:);

conf=csvread(configfile);

roi_crop=conf(3,:);
sides_cm=conf(18,1:2);

calib=calib2d(roi_crop, sides_cm);
cm2pix=@(x,y) cm2pix1(x,y, calib);

for ii =1:size(Xh,1)
    k=Xh(ii,1);
    img=rgb2gray(imread([floc, flist(k).name]));
    xcur=Xh(ii,3);
    ycur=Xh(ii,4);
    [xcur, ycur]=cm2pix(xcur, ycur);
    img=highlightpt(img, xcur, ycur);
    imwrite(img, [floc, 'trk_', flist(k).name], 'jpeg');
end



function img=highlightpt(img, xcur, ycur)

[h, w]=size(img(:,:,1));


% set the circular region for each mosquito
npts=200; th=linspace(0,2*pi, npts);
hl_rad=0:ceil(w/30);
[HLR TH]=meshgrid(hl_rad, th);
regionx0=ceil(HLR.*cos(TH));
regiony0=ceil(HLR.*sin(TH));

for jj=1:numel(xcur)
    region_ctr=[xcur(jj), ycur(jj)]';
    regionx=ceil(region_ctr(1))+regionx0;
    regiony=ceil(region_ctr(2))+regiony0;
    hlght=unique([regionx(:), regiony(:)], 'rows');

    if region_ctr(1)==0, keyboard; end

    hlght(:,2)=min(h,hlght(:,2));
    hlght(:,2)=max(1,hlght(:,2));
    hlght(:,1)=min(w,hlght(:,1));
    hlght(:,1)=max(1,hlght(:,1));
    % brighten up this part in the image
    imgind=sub2ind([h,w], hlght(:,2), hlght(:,1));
    [num thresh]=hist(double(img(imgind)),5);
    thresh=mean(thresh);
    img(imgind)=img(imgind)*1.1;
    for trail_data=1:numel(imgind)
        if img(imgind(trail_data))<thresh(1)
            img(imgind(trail_data))=img(imgind(trail_data))*.99;
        end
    end
end


function [u v]=cm2pix1(x,y, calib)

u= x/calib.pix2cm(1)+calib.center(1);
v= y/calib.pix2cm(2)+calib.center(2);



function calib=calib2d(roi, sides_cm)
%function [center, pix2cm]=calib2d(roi, sides_cm)
%
% roi
% sides_cm is 1x2 in cm of sides

if sides_cm(1)
    center=[roi(1)+roi(3)/2;
            roi(2)+roi(4)/2];
    

    tanksides=[roi(3), roi(4)];
    pix2cm=sides_cm./tanksides(1:2);
else
    pix2cm=ones(1,2);
    center=zeros(2,1);
end

calib.pix2cm=pix2cm;
calib.center=center;

%%%%%%%%%%%%%%%%%%%%% old code for reformatting old data
% if sides_cm(1)
%     cr=[[roi(1); roi(2)], [roi(1)+roi(3); roi(2)], ...
%         [roi(1); roi(2)+roi(4)], [roi(1)+roi(3); roi(2)+roi(4)]];
% 
%     cr(:,5)=cr(:,1);
% 
%     tanksides=sqrt(sum(diff(cr,1,2).^2,1));
%     center=[cr(1,2)+cr(1,1);
%                 cr(2,3)+cr(2,2);]/2;
% 
%     pix2cm=sides_cm./tanksides(1:2);
% else
%     pix2cm=ones(1,2);
%     center=zeros(2,1);
% end
% 
% calib.pix2cm=pix2cm;
% calib.center=center;