function img=quick_verify(k, X,  cm2pix, img)
% 

if nargin ==4 % image is provided for changing and passing
    xcur=X(X(:,1)==k,3);
    ycur=X(X(:,1)==k,4);
    if size(img, 3) >1
        img=rgb2gray(img);
    end
    [xcur, ycur]=cm2pix(xcur, ycur);
    img=highlightpt(img, xcur, ycur);
else
    nsee=10; nfut=3;
    trail=X(:,1)>max(1,k-nsee) & X(:,1)<=k+nfut;
    Xtr=X(trail, :);
    curr=Xtr(:,1)==k;

    ids=unique(Xtr(:,2));
    u=nan(nsee+nfut, numel(ids));
    v=u;
    u1=nan(1,numel(ids));
    v1=u1;
    jj=1;
    % p_fut=[]; p_rep=[];
    for ii=ids'
        iidx=Xtr(:,2)==ii;
        x=Xtr(iidx, 3);
        y=Xtr(iidx, 4);
        
        [val tidx]=sort(Xtr(iidx,1));
        x=x(tidx); y=y(tidx);
        

        [u(1:numel(x),jj) v(1:numel(x),jj)]=cm2pix(x,y);


        x1=Xtr(curr & Xtr(:,2)==ii,3);
        y1=Xtr(curr & Xtr(:,2)==ii,4);
        if ~isempty(x1)
            x1=x1(1); y1=y1(1);
            [u1(1,jj) v1(1,jj)]=cm2pix(x1,y1);
        end

        jj=jj+1;
    end
    
    xcur=Xtr(curr,3);
    ycur=Xtr(curr,4);
    [xcur, ycur]=cm2pix(xcur, ycur);
    good=find(~isnan(sum(u,1)));
    plot(u1, v1, 'bo', 'linewidth', 1.5, 'markersize', 10);
    if ~isempty(good)
        plot(u1(good), v1(good), 'go', 'linewidth', 1.5, 'markersize', 10);
    end
    plot(u,v);

    % quiver if needed
    if size(X,2)>10
        fl=ceil(sqrt(Xtr(curr,7)/3)*2);
        if ~isempty(fl), fl=fl(1); else fl=15; end
        hd=Xtr(curr,11:12)'; % 11:12 because the first two columns are time and id
        quiver(xcur, ycur, hd(1,:)'*10, hd(2,:)'*10, 0, 'r');
        shape=Xtr(curr,9:10)'; % ** 9:10 because the first two columns are time and id
        for ff=1:size(hd,2)
            wTb=[hd(:,ff), [-hd(2,ff), hd(1,ff)]', [xcur(ff), ycur(ff)]'; 0 0 1];
            xx=(-fl:fl)'; % based on 1:3 ratio between width and length
            A=[xx.^2, xx];
            yy=A*shape(:,ff);
            pixw=tra2b([xx,yy]', wTb);
            plot(pixw(1,:), pixw(2,:), 'g');
        end
    end
end

function img=highlightpt(img, xcur, ycur)

[h, w]=size(img(:,:,1));


% set the circular region for each mosquito
npts=200; th=linspace(0,2*pi, npts);
hl_rad=0:ceil(w/30);
[HLR TH]=meshgrid(hl_rad, th);
regionx0=ceil(HLR.*cos(TH));
regiony0=ceil(HLR.*sin(TH));

for jj=1:numel(xcur)
    region_ctr=[xcur(jj), ycur(jj)]';
    regionx=ceil(region_ctr(1))+regionx0;
    regiony=ceil(region_ctr(2))+regiony0;
    hlght=unique([regionx(:), regiony(:)], 'rows');

    if region_ctr(1)==0, keyboard; end

    hlght(:,2)=min(h,hlght(:,2));
    hlght(:,2)=max(1,hlght(:,2));
    hlght(:,1)=min(w,hlght(:,1));
    hlght(:,1)=max(1,hlght(:,1));
    % brighten up this part in the image
    imgind=sub2ind([h,w], hlght(:,2), hlght(:,1));
    [num thresh]=hist(double(img(imgind)),5);
    thresh=mean(thresh);
    img(imgind)=img(imgind)*1.1;
    for trail_data=1:numel(imgind)
        if img(imgind(trail_data))<thresh(1)
            img(imgind(trail_data))=img(imgind(trail_data))*.99;
        end
    end
end
