function stats_1way(file, cselect, ph)

fid = fopen(file, 'r');
tline = fgetl(fid);

% Split header
A(1,:) = regexp(tline, '\,', 'split');

% Parse and read rest of file
ctr = 1;
while(~feof(fid))
if ischar(tline) 
ctr = ctr + 1;
tline = fgetl(fid); 
A(ctr,:) = regexp(tline, '\,', 'split'); 
else
break; 
end
end
fclose(fid);

% [bin_centers, units, jnk1, jnk2, id, id2]= get_options(type);

if nargin<2, cselect=1:numel(id); ph=0; end

% anal='Variation in '; 
anal='Average ';

fprintf('Comparing ....\n');
for ii=cselect
    fprintf('%s ', id2{ii});
end
fprintf('\n');

%%%%% STATS
% ad_test=csvread(['../data/data1_', type, '.csv']);
ad_test=str2double(A(2:end,:));

[a1 b1 st]=anova1(ad_test(:,cselect),id2(cselect),'off');
if ph
    c=multcompare(st, 'alpha', .05, 'display', 'off');
    for ii=cselect
        for jj=cselect
            ii1=find(cselect==ii);
            jj1=find(cselect==jj);
            row=find(c(:,1)==ii1 & c(:,2)==jj1);

            if ~(c(row,5) > 0 & c(row,3) < 0)

                p(ii,jj)=1;
            else
                p(ii,jj)=0;
            end
        end
    end
end

if ~isempty(strfind(type, 'Distance to robot'))
    p(1,:)=0;
    ad_test(:,1)=0;
end
    
N=size(ad_test,1);

figure(2); gcf; clf;
set(gcf, 'units', 'normalized');
bar(nanmean(ad_test(:,cselect),1), .5, 'k');
hold on;
errorbar(nanmean(ad_test(:,cselect),1), nanstd(ad_test(:,cselect),[],1)./sqrt(sum(~isnan(ad_test(:,cselect)),1)), 'k.', 'linewidth', 2);
kk=1;
if ph
    signi=[];
    for ii=cselect
        for jj=ii+1:cselect(end)
            if p(ii,jj)
                signi(kk,:)=[find(cselect==ii),find(cselect==jj)];
                kk=kk+1;
            end
        end
    end
end

%%%% PLOT
if ph
if ~isempty(signi)
    hs=sigstar(num2cell(signi,2));
    set(hs, 'color', 'b');
end
end
set(gca, 'fontsize', 24);
set(gca, 'xticklabel', id2(cselect)); 
set(gca, 'ylim', [bin_centers(1), bin_centers(end)]);
ylabel([anal, type, units]);
set(gcf, 'position', [0.0599    0.1418    0.4443    0.7902]);
set(gca, 'position', [0.1114    0.0646    0.8710    0.9175]);
box off;
text(1, bin_centers(end)*.9, sprintf('p=%.3f', a1), 'fontsize', 25);
% title(sprintf('one-way anova: %.5f', a1));
fprintf('ANOVA=%.5f \t F=%.3f\n', a1, b1{2,5});

% resp=input('Save? []=no, other=yes: ');
% if ~isempty(resp)
%     saveas(gcf, ['/Users/sachit/Dropbox/Research/bbp_2013b/figures/', acr, '_all.eps'], 'eps');
%     saveas(gcf, ['/Users/sachit/Documents/Publish/2013/papers/bbp_2013b/supplementary/', acr, '_all.fig'], 'fig');
% end

