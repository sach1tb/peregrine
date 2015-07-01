function cleanup_data(csvfile)
% function cleanup_data(csvfile)

fprintf('keep a backup copy of the file!!\n');

csvdata=csvread(csvfile);

numFrames=max(csvdata(:,1));

% %%% shave off last five minutes
% ff=csvdata(csvdata(:,2)==1 & csvdata(:,1)>0,1);
% ff=min(ff);
% nfr=7200;
% 
% csvdata=csvdata(csvdata(:,1)>=ff & csvdata(:,1)<=ff+nfr,:);

% remove zero ids
ids=unique(csvdata(:,2));
ids2=[];
ids=ids(ids~=0);

% remove tracks that are less than 1% of the number of frames
fprintf('removing tracks less than %d frames long..\n', ...
            floor(numFrames*.01));
numIds=numel(ids);
for jj=ids'
    if sum(csvdata(:,2)==jj,1) > floor(numFrames*.01)
        ids2=[ids2, jj];
    end
end
fprintf('tracks removed=%d\n', numIds-numel(ids2));
ids=ids2;
ids2=[];

fprintf('removing tracks that do not move at all...\n')
resp=input('confirm? []=no, other=yes:');
if ~isempty(resp)
    numIds=numel(ids);
    for jj=ids
        velocity=csvdata(csvdata(:,2)==jj, 5:6);
        speed=mean(sum(velocity.^2,1));
        if speed>0.01 % cm/s
            ids2=[ids2, jj];
        end
    end
    fprintf('tracks removed=%d\n', numIds-numel(ids2));
    ids=ids2;
    ids2=[];
end


fprintf('keeping tracks longer than %d frames ...\n', ...
            floor(numFrames*.75));
resp=input('confirm? []=no, other=yes:');
if ~isempty(resp)
    numIds=numel(ids);
    for jj=ids
        if sum(csvdata(:,2)==jj,1) > floor(numFrames*.75)
            ids2=[ids2, jj];
        end
    end
    fprintf('tracks removed=%d\n', numIds-numel(ids2))
    ids=ids2;
    ids2=[];
end


if ~isempty(ids)
    csvdata=csvdata(ismember(csvdata(:,2), ids),:);
    [~, idx]=sort(csvdata(:,1));
    csvdata=csvdata(idx,:);
    if ~isempty(csvdata)
        csvwrite(csvfile, csvdata);
    else
        delete(csvfile);
    end
else
    fprintf('[!] None of the ids have enough frames\n');
end