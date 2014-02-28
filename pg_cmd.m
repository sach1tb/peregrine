function pg_cmd(frmloc, id, live, copy_config)
%function pg_cmd(frmloc, id, live)
% command line version of peregrine. Use it after setting the thresholds to
% redo the tracking, or to track for setups that do not change after
% copying the preferences...
%
% frmloc is the location of frames
% id is the image id without numbers
% live is 1 for showing images
% copy_config is the file that must be copied
%
% Sachit Butail, NYU, 4/2013 
%

global live1
global k
global roi_crop
global roi_cut
global view_type


live1=live;
view_type=3;
bgyes=0;

frmloc=[frmloc, '/'];

offline=1;
flist=dir([frmloc, id, '*.*']);
filename=flist(3).name;


[conf, vid, getfrm, ...
 nfrm fps, nt, ...
 bg, alpha, smooth blur, bb_blur, bb_size, ...
 suff max_val imh imw fgislight trktype, ...
 record_verify roi_crop roi_cut img_t, area_t, ...
 bin_t sides_cm circ, shape, split]= init_setup(frmloc, offline, filename);

if nargin ==4
    config_file=[frmloc, 'config_', suff, '.csv'];
    copyfile(copy_config, config_file);
    
    [conf, vid, getfrm, ...
     nfrm fps, nt, ...
     bg, alpha, smooth blur, bb_blur, bb_size, ...
     suff max_val imh imw fgislight trktype, ...
     record_verify roi_crop roi_cut img_t, area_t, ...
     bin_t sides_cm circ, shape, split]= init_setup(frmloc, offline, filename);
end

write_tracks=1;

%%%% file for saving tracks
csvfile=[frmloc, 'X00_', suff, '.csv'];
if exist(csvfile, 'file')
    error('[!] datafile alreadly present. the new file will overwrite it! Delete it manually first');
else
    datfile=zeros(nfrm, 20);
end

sz1=zeros(1,1000);
X=[]; P=[];

frame(nfrm,1).Zk=[];
frame(nfrm,1).Mk=[];
frame(nfrm,1).X=[];
frame(nfrm,1).P=[];
tic
for k=1:nfrm
    
    if live1   
        figure(1); gcf; clf;    
    end

    [bg nZ X P datfile cm2pix sz1 mesg frame]=onestep(sides_cm, img_t, bin_t, area_t, ...
                                        fgislight, circ, smooth, bg, alpha,...
                                        blur, bb_blur, bb_size, X, P, fps, ...
                                        datfile, getfrm, bgyes, trktype, ...
                                        write_tracks, record_verify, split, shape, sz1, frame);
end
toc

if sum(datfile(:,1))
    csvwrite([frmloc, 'X00_', suff, '.csv'], datfile);
end