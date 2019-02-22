function conf =read_conf(cfile)


% change file content
fid=fopen(cfile);
C=textscan(fid, '%s %s %f %f %f %f', 'Delimiter', ',');
for ii=1:size(C{1},1)
    fieldname=C{2}{ii};
    for jj=1:4
        if ~isnan(C{jj}(1))
            fieldvalue(jj)=C{jj}(1);
        end
    end
    conf.(fieldname)=fieldvalue;
end
for field=fieldnames(conf)'
    fprintf(fid, ', %s,', field{1});
    for jj=1:numel(conf.(field{1}))
        fprintf(fid, '%f,', conf.(field{1}));
    end
    for kk=jj:4
        fprintf(fid, '%s,', 'NaN');
    end
    fprintf(fid, '\n');
end
