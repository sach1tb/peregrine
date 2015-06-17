function [conf, vid, getfrm, nfrm, fps, nt, bg, alpha, ...
    smooth, blur, bb_blur, bb_size, suff, max_val,...
    imh, imw, fgislight, trktype, record_verify, ...
    roi_crop, roi_cut, img_t, area_t, bin_t, sides_cm, circ, shape,split]=init_setup(frmloc, offline, filename)

switch offline
    case 1 % images
        fps=0;
        if ~isempty(filename)
            ln=numel(filename);
            id=[filename(1), '*.*', filename(ln-2:ln)];
        else
            id='raw*.*'; % if you pass a direct argument
        end
        frmlist = dir([frmloc, id]);
        
        try
            cat(1,frmlist.name);
        catch me_err
            fprintf('[!] Improper filenames... all filenames to be of the same length.\n');
            fprintf('[I] Run padzero. Use help padzero to see how to run first\n');
        end
        
        nfrm=size(frmlist,1);
        getfrm=@(k) imread([frmloc, frmlist(k).name]);
        img=getfrm(1);
        suff=frmlist(1).name(1:end-6);
        vid=[];
    case 2 % video
        fps=0;
        try
        obj=VideoReader([frmloc, filename]);
        catch
%             close(gcf);
            error(sprintf('[!] Unable to read video.\n Try running on frames extracted from the video.\n*****\n'));
        end
        getfrm=@(k) read(obj, k);
        nfrm=obj.NumberOfFrames;
        img=getfrm(1);
        suff=filename(1:end-4);
        vid=[];
    case 3 % webcamera
        fps=30;
        vp=videoParams;
        vid=videoinput(vp.adaptername, vp.deviceid, vp.format);
        set(vid,'TriggerRepeat',Inf);
        vid.FrameGrabInterval =30*1/fps; % 30fps * dt
        vid.ReturnedColorspace = 'rgb';
        getfrm=@(k)  liveread(k, vid, 0);
        start(vid);
        img=getfrm(1);
        stop(vid);
        suff=datestr(now,'yyyy-mm-ddTHHMM');
        nfrm=0;
end

[imh, imw, imd]=size(img);

if exist([frmloc, 'config_', suff, '.csv'], 'file')
    conf=csvread([frmloc, 'config_', suff, '.csv']);
else
    conf=defconfig;
    conf(3,:)=[1 1 imw, imh];
end
if ~fps, fps=conf(1,3); end
if ~fps, error('check configuration'); end
nt=conf(1,1);

img_t=conf(2,1);
area_t=conf(2,2);
bin_t=conf(2,3);

roi_crop=conf(3,:);
roi_cut=conf(4:13,:);  

smooth=conf(14,:);
blur=conf(15,:);
bb_blur=conf(16,:);

fgislight=conf(17,1);
trktype=conf(17,2);
bb_size=conf(17,3);
record_verify=conf(17,4);

sides_cm=conf(18,1:2);
split=conf(18,3);
circ=conf(18,4);

shape=conf(19,1);
alpha=conf(19,2);

%%%% 
if ~bb_size, bb_size=200; end

imtype=class(img);
max_val=2^str2double(imtype(5:end))-1;
if size(img, 3) > 1, img=rgb2gray(img); end
bglist=dir([frmloc, 'background_*.*']);
if ~isempty(bglist)
    bg=imread([frmloc, bglist(1).name]);
else
    bg=img;
end


