function cleanup_data(csvfile)
% function cleanup_data(csvfile)

fprintf('Making a backup copy of the file...\n');
copyfile(csvfile, [csvfile(1:end-4), '_', datestr(now, 'ddmmyyyy_HHMMSS'), '.csv']);

csvdata=csvread(csvfile);


% %%% shave off last five minutes
% ff=csvdata(csvdata(:,2)==1 & csvdata(:,1)>0,1);
% ff=min(ff);
% nfr=7200;
% 
% csvdata=csvdata(csvdata(:,1)>=ff & csvdata(:,1)<=ff+nfr,:);

% remove extra rows
csvdata=csvdata(csvdata(:,2)~=0,:);
numFrames=max(csvdata(:,1))-min(csvdata(:,1));


% get all the target ids
ids=unique(csvdata(:,2));
ids2=[];

% remove tracks that are less than 1% of the number of frames
fprintf('Remove tracks that less than %d frames long..\n', ...
            min(5, floor(numFrames*.01)));
resp=input('confirm? []=no, other=yes:');
if ~isempty(resp)
    numIds=numel(ids);
    for jj=ids'
        if sum(csvdata(:,2)==jj,1) > floor(numFrames*.01)
            ids2=[ids2, jj];
        end
    end
    fprintf('tracks removed=%d\n', numIds-numel(ids2));
    ids=ids2';
    ids2=[];
end

fprintf('Remove tracks that do not move at all...\n')
resp=input('confirm? []=no, other=yes:');
if ~isempty(resp)
    numIds=numel(ids);
    for jj=ids'
        pos=csvdata(csvdata(:,2)==jj, 3:4);
        dist_moved=sum(sqrt((sum((diff(pos, 1)).^2,2))));
        if dist_moved>1 % cm
            ids2=[ids2, jj];
        end
    end
    fprintf('tracks removed=%d\n', numIds-numel(ids2));
    ids=ids2';
    ids2=[];
end

fprintf('Keep only tracks that are longer than %d frames ...\n', ...
            floor(numFrames*.75));
resp=input('confirm? []=no, other=yes:');
if ~isempty(resp)
    numIds=numel(ids);
    for jj=ids'
        if sum(csvdata(:,2)==jj,1) > floor(numFrames*.75)
            ids2=[ids2, jj];
        end
    end
    fprintf('tracks removed=%d\n', numIds-numel(ids2))
    ids=ids2';
    ids2=[];
end

fprintf('Remove tracks which belong to targets that are too big or small...\n')
resp=input('confirm? []=no, other=yes:');
if ~isempty(resp)
    numIds=numel(ids);
    avg_size=mean(csvdata(:,7));
    std_size=std(csvdata(:,7));
    for jj=ids'
        size_id=mean(csvdata(csvdata(:,2)==jj, 7));
        if abs(size_id-avg_size)<3*std_size
            ids2=[ids2, jj];
        end
    end
    fprintf('tracks removed=%d\n', numIds-numel(ids2));
    ids=ids2';
    ids2=[];
end

if ~isempty(ids)
    csvdata=csvdata(ismember(csvdata(:,2), ids),:);
    % sort according to time
    [~, idx]=sort(csvdata(:,1));
    csvdata=csvdata(idx,:);
    if ~isempty(csvdata)
        csvwrite(csvfile, csvdata);
        fprintf('%% manual repair done: %.3f\n', ...
                sum(csvdata(:,9)==1)/size(csvdata,1)*100);
    else
        delete(csvfile);
    end
else
    fprintf('[!] None of the ids have enough frames\n');
end