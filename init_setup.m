function [conf, vid, getfrm, nfrm, bg, suff, errmesg, maxIntensity]=init_setup(frmloc, offline, filename)
errmesg='';

switch offline
    case 1 % images
        fps=0;
        if ~isempty(filename)
            ln=numel(filename);
            % this will store the extension as well
            % assumes that the extension is a 3 letters
            id=[filename(1), '*.', filename(ln-2:ln)];
        else
            id='raw*.*'; % if you pass a direct argument
        end
        frmlist = dir([frmloc, id]);
        
        try
            cat(1,frmlist.name);
        catch me_err
            errmesg=sprintf('[!] Improper filenames... all filenames to be of the same length.\nClose application and remove frames with different file names or run padzero.\n');
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

% upgrade to new csv file style
oldfile=[frmloc, 'config_', suff, '.csv'];
newfile=[frmloc, 'config_', suff, '.mat'];
if exist(oldfile, 'file')
    conf=csvread(oldfile);
    conf=update_conf(conf, oldfile);
elseif exist(newfile, 'file')
    load(newfile, 'conf');
    
    %!-- 7/2017 adding t_gate and data association as variables 
    % that can be modified
    if ~isfield(conf, 't_gate')
        conf.t_gate=16;
    end

    if ~isfield(conf, 'da_type')
        conf.da_type=1;
    end
    %--!
    
else
    smooth=struct('hsize', [3 3], 'sigma', 0.5);
    conf=struct('nt', 1, 'fps', 30, 'img_t', 80, 'area_t', 0, 'bin_t', 0, ...
                'roi_crop', [1, 1, imw, imh], 'roi_cut', zeros(10,4), 'smooth', smooth, 'blur', zeros(1,4), 'bb_blur', zeros(1,4), 'fgislight', 1, ...
                'trktype', 1, 'bb_size', 0, 'record_verify', 0, 'sides_cm', 10*ones(1,2), ...
                'split', 0, 'circ', 1, 'shape', 0, 'alpha', 0, 'bgtype', 1, 't_gate', 16, 'da_type', 1);
end

% some datasets appear to have the tracker type 0
if conf.trktype==0, conf.trktype=1; end

imtype=class(img);
maxIntensity=2^str2double(imtype(5:end))-1;
if size(img, 3) > 1, img=rgb2gray(img); end
bglist=dir([frmloc, 'background_*.*']);
if ~isempty(bglist)
    bg=imread([frmloc, bglist(1).name]);
else
    bg=img;
end


function conf=update_conf(conf, oldfile)

conf1.nt=conf(1,1); 

conf1.img_t=conf(2,1);
conf1.area_t=conf(2,2);
conf1.bin_t=conf(2,3);

conf1.roi_crop=conf(3,:);
roi_cut=conf(4:13,:);
roi_cut=roi_cut(roi_cut(:,1)~=0,:);
conf1.roi_cut=roi_cut;
% for jj=1:size(roi_cut,1)
%     fldname=['roi_cut', sprintf('%.2d', jj)];
%     conf1.(fldname)=roi_cut(jj,:);
% end

conf1.smooth.hsize=conf(14,1:2);
conf1.smooth.sigma=conf(14,3);
conf1.blur=conf(15,:);
conf1.bb_blur=conf(16,:);

conf1.fgislight=conf(17,1)+1; % since we now have a list
conf1.trktype=conf(17,2);
conf1.bb_size=conf(17,3);
conf1.record_verify=conf(17,4);

conf1.sides_cm=conf(18,1:2);
conf1.split=conf(18,3);
conf1.circ=conf(18,4)+1;% since we now have a list

conf1.shape=conf(19,1);
conf1.alpha=conf(19,2);
if conf1.alpha
    conf1.bgtype=2;
else
    conf1.bgtype=1;
end
conf1.fps=conf(1,3);

if ~conf1.fps, error('check configuration'); end

%%%% 
if ~conf1.bb_size, conf1.bb_size=200; end

if ~conf1.trktype, conf1.trktype=1; end

%!-- backwards compatibility t_gate and da_type 8/2017
conf1.t_gate=16;
conf1.da_type=1;
%-!

conf=conf1;

% change file content
% fid=fopen(oldfile);
% for field=fieldnames(conf)'
%     fprintf(fid, ', %s,', field{1});
%     for jj=1:numel(conf.(field{1}))
%         fprintf(fid, '%f,', conf.(field{1}));
%     end
%     for kk=jj:4
%         fprintf(fid, '%s,', 'NaN');
%     end
%     fprintf(fid, '\n');
% end

% change file name

newfile=strrep(oldfile, '.csv', '.mat');
movefile(oldfile, newfile);
save(newfile, 'conf'); 

