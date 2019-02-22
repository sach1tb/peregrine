function recon3d

% os dependent
if ispc
    fileDelimiter='\';
elseif ismac || isunix
    fileDelimiter='/';
end

% get the camera calibration file 
[filename, pathname, ~]=uigetfile(...
            {'*.mat'},...
            'Choose the calibration file (camera_calibration.mat)');
load([pathname, fileDelimiter, filename]);

numberOfCameras=size(cams,2);

if numberOfCameras ~= 2
    error('this script is for 2 views only!');
end

% some fixes for old files
% for c_ii=1:numberOfCameras
%     cams(c_ii).km=[cams(c_ii).fc(1), cams(c_ii).alpha_c*cams(c_ii).fc(1), cams(c_ii).cc(1);
%                  0, cams(c_ii).fc(2), cams(c_ii).cc(2);
%                  0, 0, 1];       
% end


for cc=1:numberOfCameras
    
    figure(1); gcf; clf;

    [datafile, datapath, ~]=uigetfile(...
            {'*.csv'},...
            sprintf('Select the data file for camera %d', cc));
        
    configfile=[datapath, fileDelimiter, 'config_', datafile(5:end)];
    if ~exist(configfile, 'file')
        error('Configuration file should be in the same location');
    end
    conf=csvread(configfile);
    calib=calib2d(conf(3,:), conf(18,1:2));
    

    imglist=dir([datapath, fileDelimiter, datafile(5:end-4),'*.*']);
    if isempty(imglist)
        error('Image files should be in the same location');
    end
    img{cc}=imread([datapath, fileDelimiter, imglist(1).name]);
    
    % clean up the files before starting
    cleanup_data([datapath, fileDelimiter, datafile]);
    X00=csvread([datapath, fileDelimiter, datafile]);
    
    % only works for id==1, i.e., files should be cleaned
    rv=X00(X00(:,2)==1,3:4)';
    for jj=1:2
        rv(jj,:)=sma(rv(jj,:), 5);
    end

    % show images. for each view, extract the corner points of roi, also
    % extract the pixel positions of the trajectory back from cm using the
    % 2d calibration routine
    imshow(img{cc});
    hold on;
    
    title('Draw a rectangle on the outer surface of ROI', ...
            'fontsize', 16);
    outer=getrect(1);
    rectangle('Position', outer, 'LineStyle', '--', 'EdgeColor', 'red')
    title('Draw a rectangle on the inner surface of ROI', ...
        'fontsize', 16);
    inner=getrect(1);
    if cc==1
        x=[outer(1); outer(1)+outer(3); outer(1)+outer(3); outer(1); ...
            inner(1); inner(1)+inner(3); inner(1)+inner(3); inner(1)];
        y=[outer(2); outer(2); outer(2)+outer(4); outer(2)+outer(4);...
           inner(2); inner(2); inner(2)+inner(4); inner(2)+inner(4);];
    elseif cc==2
        x=[inner(1); inner(1)+inner(3); outer(1)+outer(3); outer(1); ...
          inner(1); inner(1)+inner(3); outer(1)+outer(3); outer(1);];
        y=[inner(2); inner(2); outer(2); outer(2);...
           inner(2)+inner(4); inner(2)+inner(4); outer(2)+outer(4); outer(2)+outer(4);];
    end
    roiPix{cc}=[x,y];
    
    % get the trajectory back in pixels for each view
    [us, vs]=cm2pix1(rv(1,:), rv(2,:), calib);
    trajectoryPix{cc}=[us; vs];

    idx1{cc}=unique(X00(X00(:,2)==1,1));
end

% find the number of points that are common to both views.
idx=intersect(idx1{1}, idx1{2});
for cc=1:numberOfCameras
    trajectoryPix{cc}=trajectoryPix{cc}(:,idx);
end

nfrm=min(size(trajectoryPix{1},2), size(trajectoryPix{2},2));


%%% improve camera parameter
roiCornerPix=[roiPix{1}'; roiPix{2}'];
npts=size(roiCornerPix,2);

% for each of the eight corner find an initial guess of 3D position
for jj=1:8
    [roiCorners(:,jj), err1]=lsTriangulate([roiPix{1}(jj,:)', ...
        roiPix{2}(jj,:)'], cams);
end

% the full parameters are camera euler angles, position, and corners. note
% that intrinsic properties are assumed fine
ea=SpinCalc('DCMtoEA321', cams(2).cTw(1:3,1:3), 'tol', .0001);
ea(ea>90 | ea<-90)=ea(ea>90 | ea<-90)-360;
% 3 euler angles, 3d position, and 8 corner positions
X0=[ea'; cams(2).cTw(1:3,4); roiCorners(:)];
range=[10 10 10 50 50 50 100*ones(1,npts*3)]';
Xmin=X0-range;
Xmax=X0+range;
options = optimset('Display','off');
warning off
[X, fval]=fmincon(@(X) cost1(X,roiCornerPix,cams), X0, [],[],[],[],Xmin,Xmax,[], options); % @(X) con1(X)
warning on;
% X=X0; fval=1;
fprintf('# of frames=%d, fval=%.1f\n',nfrm, fval);
% resp=input('proceed []=yes: ');

tank=reshape(X(7:end), 3, 8);
X0(1:6)'-X(1:6)'

cams(2).cTw(1:3,1:3)=SpinCalc('EA321toDCM', X(1:3,1)', 'tol', 0.001);

cams(2).cTw(1:3,4)=X(4:6,1);


% solve for all the rest
% for jj=1:min(size(ps,2), size(pt,2));
%     m=[ps(:,jj), pt(:,jj)];
%     rmin=r0-500; rmax=r0+500;
%     [err(jj), r(:,jj)]=localize_lsq(r0, m, get_cam_calib,[1 2]', rmin, rmax);
%     r0=r(:,jj);
% end
for jj=1:nfrm
    [r(:,jj), err(jj)]=lsTriangulate([trajectoryPix{1}(:,jj), trajectoryPix{2}(:,jj)], cams);
end



for jj=1:3
    r(jj,:)=sma(r(jj,:), 25);
end


% confirm
% figure(1); gcf; clf;
% ha=tight_subplot(1, 3, [0.01 0.01], [0.01 0.01], [0.01 0.01]);
figure(1); gcf; clf;
for cc=1:numberOfCameras
    subplot(1,numberOfCameras+1,cc); gca; cla;
%     figure(cc); gcf; clf;
    reconPix=w2cam(r, cams(cc));
    imshow(img{cc});
    hold on;
    if cc==1
    title(sprintf('Camera %d. 3D projected track (blue dashed). 2D track (red)', cc));
    else
        title(sprintf('Camera %d', cc));
    end
    
    plot(trajectoryPix{cc}(1,:), trajectoryPix{cc}(2,:), 'r');
    plot(reconPix(1,:), reconPix(2,:), 'b--');
    
    plot(roiPix{cc}(:,1), roiPix{cc}(:,2), 'wo', 'markersize', 8, 'linewidth', 2);
    tpix=w2cam(tank, cams(cc));
    
    plot(tpix(1,:), tpix(2,:), 'gx', 'markersize', 8, 'linewidth', 2);
end

% figure(1); gcf; 
% resp=input('Check the plots, and press enter: ');

% subplot(1,3,3); gca; cla;
subplot(1,numberOfCameras+1,numberOfCameras+1); gca; cla;
% figure(cc+1); gcf;
plot3(r(1,:), r(2,:), r(3,:), 'k'); hold on;
% plotcube([norm(tank(:,2)-tank(:,3)) norm(tank(:,1)-tank(:,2)) norm(tank(:,4)-tank(:,8))], tank(:,8)', .1, [0 0 1]); 
% l x w x h
% plotcube([norm(tank(:,1)-tank(:,2))  norm(tank(:,4)-tank(:,8)) norm(tank(:,2)-tank(:,3)) ], tank(:,1)', .1, [0 0 1]);
mycube2(tank');
for cc=1:numberOfCameras
    cTw=cams(cc).cTw;
    drawCam([75,30], 300, inv(cTw), 1);
end
axis image
axis off
% set(gcf, 'units', 'normalized');
% set(gcf, 'position', [0.1307    0.6096    0.6135    0.2977]);
% subplot(1,3,3); gca;
% set(gca, 'cameraposition', [5.3099    3.6842    2.0752]*1e+3);
% view(105, 7);

resp=input('Verify 3D reconstruction. Save data? []=no, other=yes: ');
if ~isempty(resp) 
   id=datafile(5:end-4);
   csvwrite([datapath, fileDelimiter, 'data_3D_', id, '.csv'], ...
       [(1:size(r,2))', ones(size(r,2),1), r']);
   set(gcf, 'paperpositionmode', 'auto');
   print('-dpng', [datapath, fileDelimiter, 'check_recon_' id, '.png']);
%    save([floc, 'calib_', id, '.mat'], 'cams', 'tank', 'roiPix'); 
end


function [u, v]=cm2pix1(x,y, calib)

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


function penalty=cost1(X, u, cams)
% cams(1).km(1,1)=X(1);
% cams(1).km(2,2)=X(2);
% 
% cams(2).km(1,1)=X(3);
% cams(2).km(2,2)=X(4);


cams(2).cTw(1:3,1:3)=SpinCalc('EA321toDCM', X(1:3,1)');

cams(2).cTw(1:3,4)=X(4:6,1);


% the corners
r=reshape(X(7:end,1), [3, size(u,2)]);

penalty=0;

for jj=1:size(r,2)

    trajectoryPix=w2cam(r(:,jj), cams(1));
    pix2=w2cam(r(:,jj), cams(2));
    
    % penalty is the distance in pixels between actual measurement and
    % projected value
    penalty=penalty+sum((u(1:2,jj)-trajectoryPix(1:2)).^2 + ...
            (u(3:4,jj)-pix2(1:2)).^2);
end

function [c, ceq]=con1(X)

c=[];

tank=reshape(X(7:end), 3, 8);

% the tank should be a cube
ls1=sqrt(sum((tank(:,1)-tank(:,2)).^2));
ls2=sqrt(sum((tank(:,3)-tank(:,4)).^2));
ls3=sqrt(sum((tank(:,5)-tank(:,6)).^2));
ls4=sqrt(sum((tank(:,7)-tank(:,8)).^2));

ss1=sqrt(sum((tank(:,2)-tank(:,3)).^2));
ss2=sqrt(sum((tank(:,2)-tank(:,4)).^2));
ss3=sqrt(sum((tank(:,5)-tank(:,8)).^2));
ss4=sqrt(sum((tank(:,7)-tank(:,6)).^2));

% ceq=abs(ls1-560) + abs(ls2-560) + abs(ls3-560) + abs(ls4-560) + abs(ss1-300) + abs(ss2-ss3) + abs(ss3-300) + abs(ss2-ss4);

ceq=(ls1-ls2).^2 + (ls2-ls3).^2 + (ls3-ls4).^2 + (ss1-ss2).^2 + (ss2-ss3).^2 + (ss3-ss4).^2;


function roiCornerPix = w2cam(rW, cam_str)
%function roiCornerPix = w2cam(rW, cam_str)
%
% rW:       3 x n matrix with world coordinates
% cam_str:  camera structure with cam.cTw as the transformation matrix
% 
% roiCornerPix:      2 x n image coordinates
% ref: http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/parameters.html


% get in camera frame
cr= tra2b(rW, cam_str.cTw);

% normalized coordinates
xn = cr(1,:)./cr(3,:);
yn = cr(2,:)./cr(3,:);

% distorted coordinates
rn2 =  xn.^2+yn.^2;
xd = xn.*(1+cam_str.kc(1)*rn2 +cam_str.kc(2)*rn2.^2);
yd = yn.*(1+cam_str.kc(1)*rn2 +cam_str.kc(2)*rn2.^2);

% pixel coordinates
rp(1,:) = cam_str.km(1,1)*xd + cam_str.km(1,2)*yd + cam_str.km(1,3)*1;
rp(2,:) = cam_str.km(2,1)*xd + cam_str.km(2,2)*yd + cam_str.km(2,3)*1;

% return
roiCornerPix=[rp(1,:);
    rp(2,:)];


function rB = tra2b(rA, bTa)
%function rB = tra2b(rA, bTa)
%
% Transforms point(s) rA from frame A to frame B through
% matrix transformation bTa. rA is [2 or 3]xn. bTa is 3x3 or 4x4. 
% The output rB is of the same dimension as rA.
%
% See also TRPA2B

% - SB Jan 15, 2009
% changed Jun 12, 2010 to include planar matrices

if(size(bTa,1)-1 ~= size(rA,1))
    error('[!] incompatible matrices');
end

if(size(bTa,1)==4)
    xB = bTa(1,1).*rA(1,:)+bTa(1,2).*rA(2,:)+bTa(1,3).*rA(3,:)+bTa(1,4).*1;
    yB = bTa(2,1).*rA(1,:)+bTa(2,2).*rA(2,:)+bTa(2,3).*rA(3,:)+bTa(2,4).*1;
    zB = bTa(3,1).*rA(1,:)+bTa(3,2).*rA(2,:)+bTa(3,3).*rA(3,:)+bTa(3,4).*1;

    rB=[xB
        yB
        zB];
elseif(size(bTa,1)==3)
    xB = bTa(1,1).*rA(1,:)+bTa(1,2).*rA(2,:)+bTa(1,3).*1;
    yB = bTa(2,1).*rA(1,:)+bTa(2,2).*rA(2,:)+bTa(2,3).*1;

    rB=[xB
        yB];
end
    

function mycube2(my_vertices)
% function mycube2(vertices)
% vertices should be in the following order
% top left, then clockwise and bottom left then clockwise
%
% ref: http://smallbusiness.chron.com/graph-cube-matlab-54144.html

my_faces= [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];

patch('Vertices', my_vertices, 'Faces', my_faces, 'FaceColor', 'b', 'FaceAlpha', .1);

function [r err]= lsTriangulate(Uk, cams)
%function [r err]= lsTriangulate(Uk, cams)
%
% Uk is 2xnc matrix



% use top and side views to get pose
nc=size(Uk,2);
con=zeros(nc*2,4);


for zz=1:nc
    % least squares triangulation
    con(2*(zz-1)+1,:)= [cams(zz).km(1,1), cams(zz).km(1,2), cams(zz).km(1,3)-Uk(1,zz)]*cams(zz).cTw(1:3,1:4);
    con(2*(zz-1)+2,:)= [cams(zz).km(2,1), cams(zz).km(2,2), cams(zz).km(2,3)-Uk(2,zz)]*cams(zz).cTw(1:3,1:4);
end

hh_A= con;
b= -hh_A(:,4);
A= hh_A(:,1:3);


r=(A'*A)\A'*b;

err=zeros(2,nc);
for zz=1:nc
    reproj=cams(zz).km*cams(zz).cTw(1:3,1:4)*[r;1];
    err(:,zz)=Uk(:,zz)-reproj(1:2)/reproj(3);
end
err=norm(err);

function OUTPUT=SpinCalc(CONVERSION,INPUT,tol,ichk)
%Function for the conversion of one rotation input type to desired output.
%Supported conversion input/output types are as follows:
%   1: Q        Rotation Quaternions
%   2: EV       Euler Vector and rotation angle (degrees)
%   3: DCM      Orthogonal DCM Rotation Matrix
%   4: EA###    Euler angles (12 possible sets) (degrees)
%
%Author: John Fuller
%National Institute of Aerospace
%Hampton, VA 23666
%John.Fuller@nianet.org
%
%Version 1.3
%June 30th, 2009
%
%Version 1.3 updates
%   SpinCalc now detects when input data is too close to Euler singularity, if user is choosing
%   Euler angle output. Prohibits output if middle angle is within 0.1 degree of singularity value.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                OUTPUT=SpinCalc(CONVERSION,INPUT,tol,ichk)
%Inputs:
%CONVERSION - Single string value that dictates the type of desired
%             conversion.  The conversion strings are listed below.
%
%   'DCMtoEA###'  'DCMtoEV'    'DCMtoQ'       **for cases that involve   
%   'EA###toDCM'  'EA###toEV'  'EA###toQ'       euler angles, ### should be
%   'EVtoDCM'     'EVtoEA###'  'EVtoQ'          replaced with the proper 
%   'QtoDCM'      'QtoEA###'   'QtoEV'          order desired.  EA321 would
%   'EA###toEA###'                              be Z(yaw)-Y(pitch)-X(roll).
%
%INPUT - matrix or vector that corresponds to the first entry in the
%        CONVERSION string, formatted as follows:
%
%        DCM - 3x3xN multidimensional matrix which pre-multiplies a coordinate
%              frame column vector to calculate its coordinates in the desired 
%              new frame.
%
%        EA### - [psi,theta,phi] (Nx3) row vector list dictating to the first angle
%                rotation (psi), the second (theta), and third (phi) (DEGREES)
%
%        EV - [m1,m2,m3,MU] (Nx4) row vector list dictating the components of euler
%             rotation vector (original coordinate frame) and the Euler 
%             rotation angle about that vector (MU) (DEGREES)
%
%        Q - [q1,q2,q3,q4] (Nx4) row vector list defining quaternion of
%            rotation.  q4 = cos(MU/2) where MU is Euler rotation angle
%
%tol - tolerance value
%ichk - 0 disables warning flags
%          1 enables warning flags (near singularities)
%**NOTE: N corresponds to multiple orientations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Output:
%OUTPUT - matrix or vector corresponding to the second entry in the
%         CONVERSION input string, formatted as shown above.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Pre-processer to determine type of conversion from CONVERSION string input
%Types are numbered as follows:
%Q=1   EV=2   DCM=3   EA=4
i_type=strfind(lower(CONVERSION),'to');
length=size(CONVERSION,2);
if length>12 || length<4,   %no CONVERSION string can be shorter than 4 or longer than 12 chars
    error('Error: Invalid entry for CONVERSION input string');
end
o_type=length-i_type;
if i_type<5,
    i_type=i_type-1;
else
    i_type=i_type-2;
end
if o_type<5,
    o_type=o_type-1;
else
    o_type=o_type-2;
end
TYPES=cell(1,4);
TYPES{1,1}='Q'; TYPES{1,2}='EV'; TYPES{1,3}='DCM'; TYPES{1,4}='EA';
INPUT_TYPE=TYPES{1,i_type};
OUTPUT_TYPE=TYPES{1,o_type};
clear TYPES
%Confirm input as compared to program interpretation
if i_type~=4 && o_type~=4,  %if input/output are NOT Euler angles
    CC=[INPUT_TYPE,'to',OUTPUT_TYPE];
    if strcmpi(CONVERSION,CC)==0;
        error('Error: Invalid entry for CONVERSION input string');        
    end
else
    if i_type==4,   %if input type is Euler angles, determine the order of rotations
        EULER_order_in=str2double(CONVERSION(1,3:5));
        rot_1_in=floor(EULER_order_in/100);     %first rotation
        rot_2_in=floor((EULER_order_in-rot_1_in*100)/10);   %second rotation
        rot_3_in=(EULER_order_in-rot_1_in*100-rot_2_in*10);   %third rotation
        if rot_1_in<1 || rot_2_in<1 || rot_3_in<1 || rot_1_in>3 || rot_2_in>3 || rot_3_in>3,
            error('Error: Invalid input Euler angle order type (conversion string).');  %check that all orders are between 1 and 3            
        elseif rot_1_in==rot_2_in || rot_2_in==rot_3_in,
            error('Error: Invalid input Euler angle order type (conversion string).');  %check that no 2 consecutive orders are equal (invalid)           
        end
        %check input dimensions to be 1x3x1
        if size(INPUT,2)~=3 || size(INPUT,3)~=1
            error('Error: Input euler angle data vector is not Nx3')            
        end
        %identify singularities        
        input_size = size(INPUT); 
        N = input_size(1); 

        % Identify singularities (second Euler angle out of range) 
        EA2 = INPUT(:,2); % (Nx1) 2nd Euler angle(s) 
        ZEROS = zeros(N,1); % (Nx1) 
        ONES = ones(N,1); % (Nx1) 
        if rot_1_in==rot_3_in % Type 2 rotation (1st and 3rd rotations about same axis) 
            if any(EA2>180*ONES) || any(EA2<ZEROS)
                error('Second input Euler angle(s) outside 0 to 180 degree range') 
            elseif any(EA2>=178*ONES) || any(EA2<=2*ONES)
                if ichk==1 
                    errordlg('Warning: Second input Euler angle(s) near a singularity (0 or 180 degrees).') 
                end 
            end 
        else % Type 1 rotation (rotations about three distinct axes) 
            if any(abs(EA2)>=90*ONES) 
                error('Second input Euler angle(s) outside -90 to 90 degree range') 
            elseif any(abs(EA2)>88*ONES) 
                if ichk==1 
                    errordlg('Warning: Second input Euler angle(s) near a singularity (-90 or 90 degrees).') 
                end 
            end 
        end   
    end
    if o_type==4,   %if output type is Euler angles, determine order of rotations
        EULER_order_out=str2double(CONVERSION(1,length-2:length));
        rot_1_out=floor(EULER_order_out/100);   %first rotation
        rot_2_out=floor((EULER_order_out-rot_1_out*100)/10);    %second rotation
        rot_3_out=(EULER_order_out-rot_1_out*100-rot_2_out*10); %third rotation
        if rot_1_out<1 || rot_2_out<1 || rot_3_out<1 || rot_1_out>3 || rot_2_out>3 || rot_3_out>3,
            error('Error: Invalid output Euler angle order type (conversion string).'); %check that all orders are between 1 and 3           
        elseif rot_1_out==rot_2_out || rot_2_out==rot_3_out,
            error('Error: Invalid output Euler angle order type (conversion string).'); %check that no 2 consecutive orders are equal
        end
    end        
    if i_type==4 && o_type~=4,  %if input are euler angles but not output
        CC=['EA',num2str(EULER_order_in),'to',OUTPUT_TYPE]; %construct program conversion string for checking against user input
    elseif o_type==4 && i_type~=4,  %if output are euler angles but not input
        CC=[INPUT_TYPE,'to','EA',num2str(EULER_order_out)]; %construct program conversion string for checking against user input
    elseif i_type==4 && o_type==4,  %if both input and output are euler angles
        CC=['EA',num2str(EULER_order_in),'to','EA',num2str(EULER_order_out)];   %construct program conversion string
    end
    if strcmpi(CONVERSION,CC)==0; %check program conversion string against user input to confirm the conversion command
        error('Error: Invalid entry for CONVERSION input string');
    end
end
clear i_type o_type CC

%From the input, determine the quaternions that uniquely describe the
%rotation prescribed by that input.  The output will be calculated in the
%second portion of the code from these quaternions.
switch INPUT_TYPE
    case 'DCM'
        if size(INPUT,1)~=3 || size(INPUT,2)~=3  %check DCM dimensions
            error('Error: DCM matrix is not 3x3xN');           
        end
        N=size(INPUT,3);    %number of orientations
        %Check if matrix is indeed orthogonal
        perturbed=NaN(3,3,N);
        DCM_flag=0;
        for ii=1:N,
            perturbed(:,:,ii)=abs(INPUT(:,:,ii)*INPUT(:,:,ii)'-eye(3)); %perturbed array shows difference between DCM*DCM' and I
            if abs(det(INPUT(:,:,ii))-1)>tol, %if determinant is off by one more than tol, user is warned.
                if ichk==1,
                    DCM_flag=1;
                end
            end
            if abs(det(INPUT(:,:,ii))+1)<0.05, %if determinant is near -1, DCM is improper
                error('Error: Input DCM(s) improper');           
            end
            if DCM_flag==1,
                errordlg('Warning: Input DCM matrix determinant(s) off from 1 by more than tolerance.')
            end
        end
        DCM_flag=0;
        if ichk==1,
            for kk=1:N,
                for ii=1:3,
                    for jj=1:3,
                        if perturbed(ii,jj,kk)>tol,   %if any difference is larger than tol, user is warned.
                            DCM_flag=1;
                        end
                    end
                end
            end
            if DCM_flag==1,
                fprintf('Warning: Input DCM(s) matrix not orthogonal to precision tolerance.')
            end
        end       
        clear perturbed DCM_flag   
        Q=NaN(4,N);
        for ii=1:N,
            denom=NaN(4,1);
            denom(1)=0.5*sqrt(1+INPUT(1,1,ii)-INPUT(2,2,ii)-INPUT(3,3,ii));
            denom(2)=0.5*sqrt(1-INPUT(1,1,ii)+INPUT(2,2,ii)-INPUT(3,3,ii));
            denom(3)=0.5*sqrt(1-INPUT(1,1,ii)-INPUT(2,2,ii)+INPUT(3,3,ii));
            denom(4)=0.5*sqrt(1+INPUT(1,1,ii)+INPUT(2,2,ii)+INPUT(3,3,ii));        
            %determine which Q equations maximize denominator
            switch find(denom==max(denom),1,'first')  %determines max value of qtests to put in denominator
                case 1
                    Q(1,ii)=denom(1);
                    Q(2,ii)=(INPUT(1,2,ii)+INPUT(2,1,ii))/(4*Q(1,ii));
                    Q(3,ii)=(INPUT(1,3,ii)+INPUT(3,1,ii))/(4*Q(1,ii));
                    Q(4,ii)=(INPUT(2,3,ii)-INPUT(3,2,ii))/(4*Q(1,ii));
                case 2
                    Q(2,ii)=denom(2);
                    Q(1,ii)=(INPUT(1,2,ii)+INPUT(2,1,ii))/(4*Q(2,ii));
                    Q(3,ii)=(INPUT(2,3,ii)+INPUT(3,2,ii))/(4*Q(2,ii));
                    Q(4,ii)=(INPUT(3,1,ii)-INPUT(1,3,ii))/(4*Q(2,ii));
                case 3
                    Q(3,ii)=denom(3);
                    Q(1,ii)=(INPUT(1,3,ii)+INPUT(3,1,ii))/(4*Q(3,ii));
                    Q(2,ii)=(INPUT(2,3,ii)+INPUT(3,2,ii))/(4*Q(3,ii));
                    Q(4,ii)=(INPUT(1,2,ii)-INPUT(2,1,ii))/(4*Q(3,ii));
                case 4
                    Q(4,ii)=denom(4);
                    Q(1,ii)=(INPUT(2,3,ii)-INPUT(3,2,ii))/(4*Q(4,ii));
                    Q(2,ii)=(INPUT(3,1,ii)-INPUT(1,3,ii))/(4*Q(4,ii));
                    Q(3,ii)=(INPUT(1,2,ii)-INPUT(2,1,ii))/(4*Q(4,ii));
            end
        end
        Q=Q';
        clear denom
    case 'EV'  %Euler Vector Input Type
        if size(INPUT,2)~=4 || size(INPUT,3)~=1   %check dimensions
            error('Error: Input euler vector and rotation data matrix is not Nx4')            
        end
        N=size(INPUT,1);
        MU=INPUT(:,4)*pi/180;  %assign mu name for clarity
        if sqrt(INPUT(:,1).^2+INPUT(:,2).^2+INPUT(:,3).^2)-ones(N,1)>tol*ones(N,1),  %check that input m's constitute unit vector
            error('Input euler vector(s) components do not constitute a unit vector')            
        end
        if MU<zeros(N,1) || MU>2*pi*ones(N,1), %check if rotation about euler vector is between 0 and 360
            error('Input euler rotation angle(s) not between 0 and 360 degrees')
        end
        Q=[INPUT(:,1).*sin(MU/2),INPUT(:,2).*sin(MU/2),INPUT(:,3).*sin(MU/2),cos(MU/2)];   %quaternion
        clear m1 m2 m3 MU
    case 'EA'        
        psi=INPUT(:,1)*pi/180;  theta=INPUT(:,2)*pi/180;  phi=INPUT(:,3)*pi/180;
        N=size(INPUT,1);    %number of orientations
        %Pre-calculate cosines and sines of the half-angles for conversion.
        c1=cos(psi./2); c2=cos(theta./2); c3=cos(phi./2);
        s1=sin(psi./2); s2=sin(theta./2); s3=sin(phi./2);
        c13=cos((psi+phi)./2);  s13=sin((psi+phi)./2);
        c1_3=cos((psi-phi)./2);  s1_3=sin((psi-phi)./2);
        c3_1=cos((phi-psi)./2);  s3_1=sin((phi-psi)./2);
        if EULER_order_in==121,
            Q=[c2.*s13,s2.*c1_3,s2.*s1_3,c2.*c13];
        elseif EULER_order_in==232,
            Q=[s2.*s1_3,c2.*s13,s2.*c1_3,c2.*c13];
        elseif EULER_order_in==313;
            Q=[s2.*c1_3,s2.*s1_3,c2.*s13,c2.*c13];
        elseif EULER_order_in==131,
            Q=[c2.*s13,s2.*s3_1,s2.*c3_1,c2.*c13];
        elseif EULER_order_in==212,
            Q=[s2.*c3_1,c2.*s13,s2.*s3_1,c2.*c13];
        elseif EULER_order_in==323,
            Q=[s2.*s3_1,s2.*c3_1,c2.*s13,c2.*c13];
        elseif EULER_order_in==123,
            Q=[s1.*c2.*c3+c1.*s2.*s3,c1.*s2.*c3-s1.*c2.*s3,c1.*c2.*s3+s1.*s2.*c3,c1.*c2.*c3-s1.*s2.*s3];
        elseif EULER_order_in==231,
            Q=[c1.*c2.*s3+s1.*s2.*c3,s1.*c2.*c3+c1.*s2.*s3,c1.*s2.*c3-s1.*c2.*s3,c1.*c2.*c3-s1.*s2.*s3];
        elseif EULER_order_in==312,
            Q=[c1.*s2.*c3-s1.*c2.*s3,c1.*c2.*s3+s1.*s2.*c3,s1.*c2.*c3+c1.*s2.*s3,c1.*c2.*c3-s1.*s2.*s3];
        elseif EULER_order_in==132,
            Q=[s1.*c2.*c3-c1.*s2.*s3,c1.*c2.*s3-s1.*s2.*c3,c1.*s2.*c3+s1.*c2.*s3,c1.*c2.*c3+s1.*s2.*s3];
        elseif EULER_order_in==213,
            Q=[c1.*s2.*c3+s1.*c2.*s3,s1.*c2.*c3-c1.*s2.*s3,c1.*c2.*s3-s1.*s2.*c3,c1.*c2.*c3+s1.*s2.*s3];
        elseif EULER_order_in==321,
            Q=[c1.*c2.*s3-s1.*s2.*c3,c1.*s2.*c3+s1.*c2.*s3,s1.*c2.*c3-c1.*s2.*s3,c1.*c2.*c3+s1.*s2.*s3];
        else
            error('Error: Invalid input Euler angle order type (conversion string)');            
        end
        clear c1 s1 c2 s2 c3 s3 c13 s13 c1_3 s1_3 c3_1 s3_1 psi theta phi
    case 'Q'
        if size(INPUT,2)~=4 || size(INPUT,3)~=1
            error('Error: Input quaternion matrix is not Nx4');            
        end
        N=size(INPUT,1);    %number of orientations 
        if ichk==1,
            if abs(sqrt(INPUT(:,1).^2+INPUT(:,2).^2+INPUT(:,3).^2+INPUT(:,4).^2)-ones(N,1))>tol*ones(N,1)
                errordlg('Warning: Input quaternion norm(s) deviate(s) from unity by more than tolerance')
            end 
        end
        Q=INPUT;
end
clear INPUT INPUT_TYPE EULER_order_in

%Normalize quaternions in case of deviation from unity.  User has already
%been warned of deviation.
Qnorms=sqrt(sum(Q.*Q,2));
Q=[Q(:,1)./Qnorms,Q(:,2)./Qnorms,Q(:,3)./Qnorms,Q(:,4)./Qnorms];

switch OUTPUT_TYPE
    case 'DCM'
        Q=reshape(Q',1,4,N);
        OUTPUT=[Q(1,1,:).^2-Q(1,2,:).^2-Q(1,3,:).^2+Q(1,4,:).^2,2*(Q(1,1,:).*Q(1,2,:)+Q(1,3,:).*Q(1,4,:)),2*(Q(1,1,:).*Q(1,3,:)-Q(1,2,:).*Q(1,4,:));
                2*(Q(1,1,:).*Q(1,2,:)-Q(1,3,:).*Q(1,4,:)),-Q(1,1,:).^2+Q(1,2,:).^2-Q(1,3,:).^2+Q(1,4,:).^2,2*(Q(1,2,:).*Q(1,3,:)+Q(1,1,:).*Q(1,4,:));
                2*(Q(1,1,:).*Q(1,3,:)+Q(1,2,:).*Q(1,4,:)),2*(Q(1,2,:).*Q(1,3,:)-Q(1,1,:).*Q(1,4,:)),-Q(1,1,:).^2-Q(1,2,:).^2+Q(1,3,:).^2+Q(1,4,:).^2];
    case 'EV'
        MU=2*atan2(sqrt(sum(Q(:,1:3).*Q(:,1:3),2)),Q(:,4));
        if sin(MU/2)~=zeros(N,1),
            OUTPUT=[Q(:,1)./sin(MU/2),Q(:,2)./sin(MU/2),Q(:,3)./sin(MU/2),MU*180/pi];
        else
            OUTPUT=NaN(N,4);
            for ii=1:N,
                if sin(MU(ii,1)/2)~=0,
                    OUTPUT(ii,1:4)=[Q(ii,1)/sin(MU(ii,1)/2),Q(ii,2)/sin(MU(ii,1)/2),Q(ii,3)/sin(MU(ii,1)/2),MU(ii,1)*180/pi];
                else
                    OUTPUT(ii,1:4)=[1,0,0,MU(ii,1)*180/pi];
                end
            end
        end
    case 'Q'
        OUTPUT=Q;
    case 'EA'
        if EULER_order_out==121,
            psi=atan2((Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)),(Q(:,2).*Q(:,4)-Q(:,1).*Q(:,3)));
            theta=acos(Q(:,4).^2+Q(:,1).^2-Q(:,2).^2-Q(:,3).^2);
            phi=atan2((Q(:,1).*Q(:,2)-Q(:,3).*Q(:,4)),(Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)));
	  Euler_type=2;
        elseif EULER_order_out==232;
            psi=atan2((Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3)),(Q(:,3).*Q(:,4)-Q(:,1).*Q(:,2)));
            theta=acos(Q(:,4).^2-Q(:,1).^2+Q(:,2).^2-Q(:,3).^2);
            phi=atan2((Q(:,2).*Q(:,3)-Q(:,1).*Q(:,4)),(Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)));
	  Euler_type=2;
        elseif EULER_order_out==313;
            psi=atan2((Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)),(Q(:,1).*Q(:,4)-Q(:,2).*Q(:,3)));
            theta=acos(Q(:,4).^2-Q(:,1).^2-Q(:,2).^2+Q(:,3).^2);
            phi=atan2((Q(:,1).*Q(:,3)-Q(:,2).*Q(:,4)),(Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3)));
	  Euler_type=2;
        elseif EULER_order_out==131;
            psi=atan2((Q(:,1).*Q(:,3)-Q(:,2).*Q(:,4)),(Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)));
            theta=acos(Q(:,4).^2+Q(:,1).^2-Q(:,2).^2-Q(:,3).^2);
            phi=atan2((Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)),(Q(:,3).*Q(:,4)-Q(:,1).*Q(:,2)));
	  Euler_type=2;
        elseif EULER_order_out==212;
            psi=atan2((Q(:,1).*Q(:,2)-Q(:,3).*Q(:,4)),(Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3)));
            theta=acos(Q(:,4).^2-Q(:,1).^2+Q(:,2).^2-Q(:,3).^2);
            phi=atan2((Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)),(Q(:,1).*Q(:,4)-Q(:,2).*Q(:,3)));
	  Euler_type=2;
        elseif EULER_order_out==323;
            psi=atan2((Q(:,2).*Q(:,3)-Q(:,1).*Q(:,4)),(Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)));
            theta=acos(Q(:,4).^2-Q(:,1).^2-Q(:,2).^2+Q(:,3).^2);
            phi=atan2((Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3)),(Q(:,2).*Q(:,4)-Q(:,1).*Q(:,3)));
	  Euler_type=2;
        elseif EULER_order_out==123;
            psi=atan2(2.*(Q(:,1).*Q(:,4)-Q(:,2).*Q(:,3)),(Q(:,4).^2-Q(:,1).^2-Q(:,2).^2+Q(:,3).^2));
            theta=asin(2.*(Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)));
            phi=atan2(2.*(Q(:,3).*Q(:,4)-Q(:,1).*Q(:,2)),(Q(:,4).^2+Q(:,1).^2-Q(:,2).^2-Q(:,3).^2));
	  Euler_type=1;
        elseif EULER_order_out==231;
            psi=atan2(2.*(Q(:,2).*Q(:,4)-Q(:,1).*Q(:,3)),(Q(:,4).^2+Q(:,1).^2-Q(:,2).^2-Q(:,3).^2));
            theta=asin(2.*(Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)));
            phi=atan2(2.*(Q(:,1).*Q(:,4)-Q(:,3).*Q(:,2)),(Q(:,4).^2-Q(:,1).^2+Q(:,2).^2-Q(:,3).^2));
	  Euler_type=1;
        elseif EULER_order_out==312;
            psi=atan2(2.*(Q(:,3).*Q(:,4)-Q(:,1).*Q(:,2)),(Q(:,4).^2-Q(:,1).^2+Q(:,2).^2-Q(:,3).^2));
            theta=asin(2.*(Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3)));
            phi=atan2(2.*(Q(:,2).*Q(:,4)-Q(:,3).*Q(:,1)),(Q(:,4).^2-Q(:,1).^2-Q(:,2).^2+Q(:,3).^2));
	  Euler_type=1;
        elseif EULER_order_out==132;
            psi=atan2(2.*(Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3)),(Q(:,4).^2-Q(:,1).^2+Q(:,2).^2-Q(:,3).^2));
            theta=asin(2.*(Q(:,3).*Q(:,4)-Q(:,1).*Q(:,2)));
            phi=atan2(2.*(Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)),(Q(:,4).^2+Q(:,1).^2-Q(:,2).^2-Q(:,3).^2));
	  Euler_type=1;
        elseif EULER_order_out==213;
            psi=atan2(2.*(Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)),(Q(:,4).^2-Q(:,1).^2-Q(:,2).^2+Q(:,3).^2));
            theta=asin(2.*(Q(:,1).*Q(:,4)-Q(:,2).*Q(:,3)));
            phi=atan2(2.*(Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)),(Q(:,4).^2-Q(:,1).^2+Q(:,2).^2-Q(:,3).^2));
	  Euler_type=1;
        elseif EULER_order_out==321;
            psi=atan2(2.*(Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)),(Q(:,4).^2+Q(:,1).^2-Q(:,2).^2-Q(:,3).^2));
            theta=asin(2.*(Q(:,2).*Q(:,4)-Q(:,1).*Q(:,3)));
            phi=atan2(2.*(Q(:,1).*Q(:,4)+Q(:,3).*Q(:,2)),(Q(:,4).^2-Q(:,1).^2-Q(:,2).^2+Q(:,3).^2));
	  Euler_type=1;
        else
            error('Error: Invalid output Euler angle order type (conversion string).');           
        end
        if(isreal([psi,theta,phi]))==0,
	        error('Error: Unreal Euler output.  Input resides too close to singularity.  Please choose different output type.')
        end
        OUTPUT=mod([psi,theta,phi]*180/pi,360);  %deg
        if Euler_type==1,
	        sing_chk=find(abs(theta)*180/pi>89.9);
	        sing_chk=sort(sing_chk(sing_chk>0));
	        if size(sing_chk,1)>=1,
		        error('Error: Input rotation #%s resides too close to Type 1 Euler singularity.\nType 1 Euler singularity occurs when second angle is -90 or 90 degrees.\nPlease choose different output type.',num2str(sing_chk(1,1)));
	        end
        elseif Euler_type==2,
	        sing_chk=[find(abs(theta*180/pi)<0.1);find(abs(theta*180/pi-180)<0.1);find(abs(theta*180/pi-360))<0.1];
	        sing_chk=sort(sing_chk(sing_chk>0));
	        if size(sing_chk,1)>=1,
		        error('Error: Input rotation #%s resides too close to Type 2 Euler singularity.\nType 2 Euler singularity occurs when second angle is 0 or 180 degrees.\nPlease choose different output type.',num2str(sing_chk(1,1)));
	        end
        end
end



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




