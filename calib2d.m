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