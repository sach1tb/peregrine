function padzero(directory, prefix)
%function padzero(directory, prefix)
%
% directory is the directory where the files are kept
% prefix is the portion of the file just before the number
% assuming that the last four letters in a filename are _dot_extension
% 
% example:
% padzero('/home/Downloads/test/test1', '20110623_8f2_cam0-6_')

prefix_length=numel(prefix);
flist=dir([directory, '/', prefix, '*']);

ext=flist(1).name(end-3:end);

fprintf('Showing a sample list of name changes ....\n\n');
pause(1);

for ii=1:100:size(flist,1)
    oldname=[directory, '/', flist(ii).name];
    name_length=numel(flist(ii).name);
    ident=flist(ii).name(end-(name_length-prefix_length)+1:end-4);
    newname=[directory, '/', sprintf('%s%.6d%s', prefix, str2double(ident), ext)];
    if ~strcmp(oldname, newname)
        fprintf('%s  -> %s\n', oldname, newname);
    end
end


resp=input('Does that look okay? ...Proceed with changing file names (keep a copy!) []=no, other=yes: ');

  
if ~isempty(resp)
    for ii=1:size(flist,1)
        oldname=[directory, '/', flist(ii).name];
        name_length=numel(flist(ii).name);
        ident=flist(ii).name(end-(name_length-prefix_length)+1:end-4);
        newname=[directory, '/', sprintf('%s%.6d%s', prefix, str2double(ident), ext)];
        if ~strcmp(oldname, newname)
            movefile(oldname, newname);
        end
    end
    fprintf('[I] check filenames and rerun... \n');
end