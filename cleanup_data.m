function cleanup_data(csvfile, nfrm)

csvdata=csvread(csvfile);
if nargin==1
    nfrm=max(csvdata(:,1));
end

% %%% shave off last five minutes
% ff=csvdata(csvdata(:,2)==1 & csvdata(:,1)>0,1);
% ff=min(ff);
% nfr=7200;
% 
% csvdata=csvdata(csvdata(:,1)>=ff & csvdata(:,1)<=ff+nfr,:);

ids=unique(csvdata(:,2));
ids2=[];
ids=ids(ids~=0);
for jj=ids'
    if sum(csvdata(:,2)==jj,1) > floor(nfrm*.75)
        ids2=[ids2, jj];
    end
end

if ~isempty(ids2)
    csvdata=csvdata(ismember(csvdata(:,2), ids2),:);
    [yy idx]=sort(csvdata(:,1));
    csvdata=csvdata(idx,:);
    if ~isempty(csvdata)
        csvwrite(csvfile, csvdata);
    else
        delete(csvfile);
    end
else
    fprintf('[!] None of the ids have enough frames\n');
end