function img = liveread(k, vid, record)

img=getdata(vid,1);
if record
    imwrite(img, sprintf('/tmp/raw_%s_%.4d.bmp', savedir, ...
            datestr(now, 'HH-MM-SS'), k), 'bmp');
end
% img=rgb2gray(img);
