function [conf] =defconfig

conf=[  5, 0, 30 0; % nt, nr, fps
        80, 0, 0, 0; % img_t, area_t, bin_t
        zeros(1,4); % crop
        zeros(10,4); % cut
        0, 0, 0, 0; % 14: Smoothing use if the images are noisy (slows down)
        0, 0, 0, 0; % 15: Blurring use if single targets are marked as several due to moving parts
        15, 15, 0, 10; % 16: Blurring the biggest blob, useful when there are small targets and a big disturbing blob...
        0, 0, 0, 0; % fgislight, trktype, bb_size, record_verify
        0, 0, 0, 0; % length, breadth, split, crop_shape
        0, 0, 0, 0; % shape, background_type, available space
        0, 0, 0, 0; 
        0, 0, 0, 0; 
        0, 0, 0, 0; 
        0, 0, 0, 0; 
        0, 0, 0, 0; 
        0, 0, 0, 0; 
        0, 0, 0, 0; 
        0, 0, 0, 0; 
        0, 0, 0, 0; 
        0, 0, 0, 0; 
        0, 0, 0, 0;];

