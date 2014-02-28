function stats_2way(type, cselect)


[bin_centers, units, jnk1, jnk2, id, id2]= get_options(type);

if nargin<2, cselect=1:numel(id); end

% anal='Variation in '; 
anal='Average ';



%%%%% STATS
data=csvread(['../data/data2_', type, '.csv']);
reps=size(data,1)/numel(id);
data2=[];
for jj=cselect
    data2=[data2; data(reps*(jj-1)+1:reps*jj,:)];
end

data2(isnan(data2))=rand;



fprintf('Rows are conditions columns is time\n');
p=anova2(data2,reps);
fprintf('Time effect=%.5f\ncondition effect=%.5f\ninteraction=%.5f\n', p(1), p(2), p(3));


% try anovan
% group 1
% g1=ones(size(data2));
% for ii=1:size(data2,2)
%     g1(:,ii)=ii;
% end
% 
% 
% % group 2
% g2=ones(size(data2));
% jj=1;
% for ii=1:reps:size(data2,1)
%     g2(ii:ii+reps-1,:)=jj;
%     jj=jj+1;
% end
% 
% anovan(data2(:), {g1(:), g2(:)})