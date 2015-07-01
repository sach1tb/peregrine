function get_camera_calib
% generate camera calib function
% this script will generate a camera calibration function for a day using
% the extrinsic calibration from the images specified. It will also show a
% mean and variance of the extrinsic parameters
%
% The following functions are called from the calibration toolbox
% extract_grid, compute_extrinsic

% clear all;

addpath('./TOOLBOX_calib');

%% change these values (or user input)

calibDate=input('Date when this calibration was done(mmmddyyyy):', 's');
numberOfCameras=input('Total number of cameras used (>1):');

% location where intrinsic calibration .mat files are kept
MIRROR=0;
for cc=1:numberOfCameras
    [camera{cc}.filename, camera{cc}.pathname, filterindex]=uigetfile(...
            {'*.mat', 'Calibration files (Calib_results.mat)'},...
            sprintf('Choose the calibration file for camera %d', cc));
    resp=input('Is this camera a mirror view? []=no, other=yes: ');  
    if ~isempty(resp)
        MIRROR=cc;
    end
end

% intrinsicCameraCalibrationLocations=[ ...
%             '../sampledata/multiview/20110221_202/calib/cam1/';
%             '../sampledata/multiview/20110221_202/calib/cam2/';
%             '../sampledata/multiview/20110221_202/calib/cam3/'];


% corresponding images for extrinsic. the identifiers are used to list the
% image
% cameraIdList=['20110221_cam1'; '20110221_cam2'; '20110221_cam3';];  

% extrinsicImageLocation=[uigetdir('', 'Select the directory with extrinsic frames'), '/'];
% extrinsicImageLocation='../sampledata/multiview/20110221_202/calib/extr/';

% if there is a mirror in the setup set this variable to that camera number
if MIRROR
    fprintf('[I] This calibration has a mirror in it. Please make sure!!!!\n');
end
    

%%% ********* No change beyond this point ********** 
% processing
% numberOfCameras=size(cameraIdList, 1);
% cam_colors=colormap(lines(size(cameraIdList,1)));
cam_colors=colormap(lines(numberOfCameras));




if ispc, tmpdir='C:\Windows\Temp\'; else tmpdir='/tmp/'; end



for ii=1:numberOfCameras
       
%     matfile=strcat(intrinsicCameraCalibrationLocations(ii,:), 'Calib_Results.mat');
    load([camera{ii}.pathname, '/', camera{ii}.filename], 'fc','cc','kc','alpha_c');
%     load(matfile, 'fc','cc','kc','alpha_c');
    
      
    % compute the extrinsic for all images available for this id
    [filename, extrinsicImageLocation, ~]=uigetfile(...
            {'*.bmp; *.tif; *.png; *.jpg; *.jpeg;', 'Image files'}, ...
            sprintf('Choose an extrinsic image file for camera %d', ii));
    cameraIdList(ii,:)=filename(1:end-8);
    extrinsicImageLocation=[extrinsicImageLocation, '/'];
    
    flist=dir(strcat(extrinsicImageLocation, cameraIdList(ii,:),'*.*'));
    ni=size(flist,1);
    
    for im_ii=1:ni
        disp(strcat('processing ', flist(im_ii).name, '....'));
        I=double(imread(strcat(extrinsicImageLocation, flist(im_ii).name)));
        if size(I,3)>1,
            I = I(:,:,2);
        end;
        wintx = 5;
        winty = 5;

        % function calls from camera calibration toolbox
        [x_ext,X_ext,n_sq_x,n_sq_y,ind_orig,ind_x,ind_y] = extract_grid(I,wintx,winty,fc,cc,kc);
        [omc_ext,Tc_ext,Rc_ext,H_ext] = compute_extrinsic(x_ext,X_ext,fc,cc,kc,alpha_c);
        
        %       for mirror reflections only
        if ii==MIRROR
%             I=I';
%             I=fliplr(I);
           Rc_ext=Rc_ext*[0 1 0;
                            1 0 0;
                            0 0 1];
%             omc_ext=rodrigues(Rc_ext);
        end
        
        [x_reproj] = project_points2(X_ext,omc_ext,Tc_ext,fc,cc,kc,alpha_c);
        err_reproj = x_ext - x_reproj;
        err_std2 = std(err_reproj')';
        
        
        fprintf(1,'\n\nExtrinsic parameters:\n\n');
        fprintf(1,'Translation vector: Tc_ext = [ %3.6f \t %3.6f \t %3.6f ]\n',Tc_ext);
        fprintf(1,'Rotation vector:   omc_ext = [ %3.6f \t %3.6f \t %3.6f ]\n',omc_ext);
        fprintf(1,'Rotation matrix:    Rc_ext = [ %3.6f \t %3.6f \t %3.6f\n',Rc_ext(1,:)');
        fprintf(1,'                               %3.6f \t %3.6f \t %3.6f\n',Rc_ext(2,:)');
        fprintf(1,'                               %3.6f \t %3.6f \t %3.6f ]\n',Rc_ext(3,:)');
        fprintf(1,'Pixel error:           err = [ %3.5f \t %3.5f ]\n\n',err_std2); 

        
        cam(ii).extr_calib(:,:,im_ii)=[Rc_ext, Tc_ext;
                    0 0 0 1];
    end
    
    
    % this part is to display the results
    Basis = [X_ext(:,[ind_orig ind_x ind_orig ind_y ind_orig ])];

    VX = Basis(:,2) - Basis(:,1);
    VY = Basis(:,4) - Basis(:,1);

    nX = norm(VX);
    nY = norm(VY);

    VZ = min(nX,nY) * cross(VX/nX,VY/nY);

    Basis = [Basis VZ];

    [x_basis] = project_points2(Basis,omc_ext,Tc_ext,fc,cc,kc,alpha_c);

    dxpos = (x_basis(:,2) + x_basis(:,1))/2;
    dypos = (x_basis(:,4) + x_basis(:,3))/2;
    dzpos = (x_basis(:,6) + x_basis(:,5))/2;



    figure(2);
    image(I);
    colormap(gray(256));
    hold on;
    plot(x_ext(1,:)+1,x_ext(2,:)+1,'r+');
    plot(x_reproj(1,:)+1,x_reproj(2,:)+1,'yo');
    h = text(x_ext(1,ind_orig)-25,x_ext(2,ind_orig)-25,'O');
    set(h,'Color','g','FontSize',14);
    h2 = text(dxpos(1)+1,dxpos(2)-30,'X');
    set(h2,'Color','g','FontSize',14);
    h3 = text(dypos(1)-30,dypos(2)+1,'Y');
    set(h3,'Color','g','FontSize',14);
    h4 = text(dzpos(1)-10,dzpos(2)-20,'Z');
    set(h4,'Color','g','FontSize',14);
    plot(x_basis(1,:)+1,x_basis(2,:)+1,'g-','linewidth',2);
    title('Image points (+) and reprojected grid points (o)');
    hold off;

    print('-dpng', sprintf('%s/extr_cam_%d.png', extrinsicImageLocation, ii));


    cam(ii).fc=fc;
    cam(ii).cc=cc;
    cam(ii).kc=kc;
    cam(ii).alpha_c=alpha_c;
    
    clear fc cc kc alpha_c
end

% fn_name=[extrinsicImageLocation, 'get_cam_calib_', calibDate, '.m'];

% reset the world frame to the first camera frame
cTw1=cam(1).extr_calib(:,:,1);

for c_ii=numberOfCameras:-1:1
    for im_ii=1:ni
        cam(c_ii).extr_calib(:,:,im_ii)=cam(c_ii).extr_calib(:,:,im_ii)/(cam(1).extr_calib(:,:,im_ii));
    end
end

% find the average of all extrinsic calibrations
for c_ii=1:numberOfCameras
    if ni>1
        std(cam(c_ii).extr_calib,0,3)      
        for im_ii=1:ni
            rvec(:,im_ii)=rodrigues(cam(c_ii).extr_calib(1:3,1:3,im_ii));
        end
        
        cam(c_ii).cTw=[rodrigues(mean(rvec,2)), mean(cam(c_ii).extr_calib(1:3,4,:),3);
                                    0       0       0       1];
    else
        cam(c_ii).cTw=cam(c_ii).extr_calib;
    end
end


% ff=fopen(fn_name,'w');
% fprintf(ff,'\nfunction cam = get_cam_calib_%s(camid)', calibDate);
% 
% fprintf(ff,'\n%% ============================ ');
% fprintf(ff,'\n%% Date of calibration:%s', calibDate);
% fprintf(ff,'\n%% ============================ ');
% 
% fprintf(ff, '\n\nswitch camid');
% 
% for c_ii=1:numberOfCameras
%     fprintf(ff,'\ncase %d', c_ii);
%     fprintf(ff,'\ncam.id=\''%s\'';', cameraIdList(c_ii,:));
%     fprintf(ff,'\ncam.km=[%f \t %f \t %f;\n%f \t %f \t %f;\n%f \t %f \t %f];',...
%                 cam(c_ii).fc(1), cam(c_ii).alpha_c*cam(c_ii).fc(1), cam(c_ii).cc(1), ...
%                  0, cam(c_ii).fc(2), cam(c_ii).cc(2), ...
%                  0, 0, 1);
%              
%     fprintf(ff,'\ncam.kc1=%f;', cam(c_ii).kc(1));
%     fprintf(ff,'\ncam.kc2=%f;', cam(c_ii).kc(2));
%     
%     fprintf(ff,'\ncam.trm=[%f \t %f \t %f \t %f;', cam(c_ii).cTw(1,:));    
%     fprintf(ff,'\n%f \t %f \t %f \t %f;', cam(c_ii).cTw(2,:));    
%     fprintf(ff,'\n%f \t %f \t %f \t %f;', cam(c_ii).cTw(3,:));    
%     fprintf(ff,'\n%f \t %f \t %f \t %f];', cam(c_ii).cTw(4,:));        
% 
%     fprintf(ff,'\ncam.color=[%.3f %.3f %.3f];', cam_colors(c_ii,1), cam_colors(c_ii,2), cam_colors(c_ii,3));
% end
% fprintf(ff, '\nend');
% fclose(ff);
    

figure(4); gcf; clf;
for c_ii=1:numberOfCameras
    % km matrix
    cam(c_ii).km=[cam(c_ii).fc(1), cam(c_ii).alpha_c*cam(c_ii).fc(1), cam(c_ii).cc(1);
                 0, cam(c_ii).fc(2), cam(c_ii).cc(2);
                 0, 0, 1];
    cam(c_ii).trm=cam(c_ii).cTw;
    cam(c_ii).id=cameraIdList(c_ii,:);
    cam(c_ii).wTc=inv(cam(c_ii).cTw);
    cam(c_ii).kc1=cam(c_ii).kc(1);
    cam(c_ii).kc2=cam(c_ii).kc(2);
    cTw=cam(c_ii).cTw;
    wTc=inv(cTw);
    drawCam([75,30], 300, cTw, 1);
    text(wTc(1,4)+15, wTc(2,4)+15, wTc(3,4)+15, ...
        sprintf('%d', c_ii), 'FontSize', 15, 'Color', cam_colors(c_ii,:));
    xlabel('x(mm)'); ylabel('y(mm)'); zlabel('z(mm)');
    hold on;
end
plotFrame(cTw1(1:3,4), cTw1(1:3,1:3), 50);
xpt=cTw1(1:3,4)+cTw(1:3,1)*50 - cTw(1:3,2)*5;
text(xpt(1), xpt(2), xpt(3), 'X', 'FontSize', 15, 'Color', 'r');
% plotting the checkerboard too. the four corners are roughly 20cm away 
cb_side=200;
checkerboard=[cTw1(1:3,4) cTw1(1:3,4)+cTw1(1:3,1)*cb_side, ...
              cTw1(1:3,4)+cTw1(1:3,1)*cb_side + cTw1(1:3,2)*cb_side, ...
              cTw1(1:3,4)+cTw1(1:3,2)*cb_side];
patch(checkerboard(1,:), checkerboard(2,:), checkerboard(3,:), [.75 .75 .75]);
view(3);
axis image
print('-dpng', sprintf('%s/cam_arr1.png', extrinsicImageLocation));
view(2);
print('-dpng', sprintf('%s/cam_arr2.png', extrinsicImageLocation));
saveas(gcf, sprintf('%s/cam_arr.fig', extrinsicImageLocation), 'fig');

fprintf('[I] Calibration file saved in %s\n Copy to calib after verifying the camera arrangement\n', extrinsicImageLocation);

cams=cam;
save([extrinsicImageLocation, '/camera_calibration.mat'], 'cams');

function drawCam(fov, fr, cTw, frame_yn)
%function drawCam(fov, fr, cTw, frame_yn)
%
% fov is [2x1] field of view matrix fov(1) is angle of view wide in
% degrees, fov(2) is angle of view vertical in degrees
% fr is range of field of view. 
% cTw is the camera transformation matrix in the world frame. cTw is either
% [ 4 x 4 ] or [3 x 3]
% frame_yn is 0 if you don't want the orthogonal frame on the camera and 1 if you
% want
%
% Example:
% drawCam([30,45], 100, cTw, 1);


gca; hold on;

% if size(cTw,1) ==3  
%     cTw1(1:2,1:2) = cTw;
%     cTw1(1:3,3)= 
% exit

% camera cylinder radius
cylr=fr/25;
[camx camy camz] = cylinder(cylr);

% camera cylinder height
camz=camz*fr/10;

aovw = fov(1)*pi/180; % angle of view wide
aovv = fov(2)*pi/180;

% camera base side
side_b=fr/15;

xbase = [ -side_b, side_b, side_b, -side_b] ;
ybase = [ -side_b, -side_b, side_b, side_b];
zbase = [ 0 , 0 , 0, 0];

% creating field of view 
my=tan(aovv/2)*fr; mx=tan(aovw/2)*fr;
xfov=[0 mx mx; 0 -mx mx; 0 -mx -mx;0 mx -mx]';
yfov=[0 my -my;0 my my; 0 my -my; 0 -my -my]';
zfov=[0 fr fr; 0 fr fr; 0 fr fr; 0 fr fr]';

xyz=cTw\eye(4);
or=cTw\[0 0 0 1]';
if (frame_yn), plotFrame(or, xyz, fr/2, 3, eye(3),0); end
[camxg, camyg, camzg]=trpa2b(camx, camy, camz, inv(cTw));
surf(camxg, camyg, camzg,'FaceColor', 'k', ...
            'FaceAlpha', 1, 'EdgeColor', 'none'); hold on;
[xfovg yfovg zfovg]=trpa2b(xfov, yfov, zfov, inv(cTw));
patch(xfovg, yfovg, zfovg, [.5 .5 .5], 'EdgeColor', ...
            [.5 .5 .5], 'FaceAlpha',0.1);

% uncomment to create base for the camera        
% [xbaseg, ybaseg, zbaseg]=trpa2b(xbase, ybase, zbase, inv(cTw));
% patch(xbaseg, ybaseg, zbaseg, [.5 .5 .5], 'EdgeColor', ...
%             [.5 .5 .5], 'FaceAlpha',0.5);
function plotFrame(p, R, s, varargin)
%plotFrame(p, R, s, varargin)
%
% p is the position of the origin
% R is the rotation matrix
% s is the scale
%
% varargin{1} is lw, linewidth
% varargin{2} is color of each axis, row of [3x3] matrix 
% varargin{3} is whether arrow (1) should be used or quiver (0)

gca; hold on;

if(nargin>3)
    lw=varargin{1};
    color=varargin{2};
else
    lw=1.5;
    color=[ 1 0 0
            0 1 0
            0 0 0];
end

if nargin==6
    if varargin{3}
        arrow3(p(1:3)*ones(1,3), p(1:3)*ones(1,3)+R(1:3,1:3)*20, 'd-1', 1, 2);
%         arrow3([0,0,0; 0,0,0; 0,0,0], [200,0,0;0,200,0;0,0,500], 'd-1', 1, 2);
    end
else
    % x 
    quiver3(p(1), p(2), p(3), R(1,1), R(2,1), R(3,1), s, 'Color', color(1,:), 'LineWidth', lw);
    % y
    quiver3(p(1), p(2), p(3), R(1,2), R(2,2), R(3,2), s, 'Color', color(2,:), 'LineWidth', lw);
    % z
    quiver3(p(1), p(2), p(3), R(1,3), R(2,3), R(3,3), s, 'Color', color(3,:), 'LineWidth', lw);
    
end
function [bx by bz] = trpa2b(ax, ay, az, bTa)
%function [bx by bz] = trpa2b(ax, ay, az, bTa)
%
% ax, ay, az can be any dimensions. 
% The output bx, by, bz are of the same dimensions. 
% bTa is a [4x4] transformation matrix
%
% See also TRA2B


bx = bTa(1,1).*ax+bTa(1,2).*ay+bTa(1,3).*az+bTa(1,4).*1;
by = bTa(2,1).*ax+bTa(2,2).*ay+bTa(2,3).*az+bTa(2,4).*1;
bz = bTa(3,1).*ax+bTa(3,2).*ay+bTa(3,3).*az+bTa(3,4).*1;
