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
    