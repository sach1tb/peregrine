function vp = videoParams

if ismac % Hold up an image or a screenshot 
    vp.adaptername='macvideo';
    vp.deviceid=1;
    vp.format='YCbCr422_1280x720';
elseif isunix
    vp.adaptername='linuxvideo';
    vp.deviceid=1;
    vp.format='YUYV_640x480';
end