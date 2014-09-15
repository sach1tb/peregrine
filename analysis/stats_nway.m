function stats_nway(type, gxg)

close all;

[bin_centers, units, jnk1, jnk2, id, id2]= get_options(type);

%%%%% STATS
data=csvread(['../data/data1_', type, '.csv']);


switch gxg
    case 'ETHxSIZ'
        cselect=[1:8];
        g1=[1 2 3 4 1 2 3 4];
        g2=[1 1 1 1 2 2 2 2];
        varnames={'Ethanol conc', 'Size shoal'};
%     case 'SPDxCFG'
%         cselect=[4 5 8 9 11 12];
%         g1=[1 1 2 2 3 3];
%         g2=[1 2 1 2 1 2];
%         varnames={'Robot Speed', 'Shoal Configuration'};    
end


fprintf('Comparing ....\n');
for ii=cselect
    fprintf('%s ', id2{ii});
end
fprintf('\n');

data2=data(:, cselect);



g1=ones(size(data2,1),1)*g1;
g2=ones(size(data2,1),1)*g2;
figure(1); gcf; clf;
set(gcf, 'units', 'normalized');
bar(nanmean(data(:,cselect),1), .5, 'k');
hold on;
errorbar(nanmean(data(:,cselect),1), nanstd(data(:,cselect),[],1)./sqrt(sum(~isnan(data(:,cselect)),1)), 'k.', 'linewidth', 2);
set(gca, 'fontsize', 10);
set(gca, 'xticklabel', id2(cselect)); 
set(gca, 'ylim', [bin_centers(1), bin_centers(end)]);
ylabel(['Average ', type, units]);
box off;

anovan(data2(:), {g1(:), g2(:)}, 'model', 'interaction', 'varnames', varnames)
resp=input('Save? []=yes, other=no: ');
if isempty(resp)
    saveas(gcf, ['../plots/', type, '_', gxg, '_nway.png'], 'png');
%     saveas(gcf, ['/Users/sachit/Documents/Publish/2013/papers/bbp_2013b/supplementary/', acr, '_all.fig'], 'fig');
end


