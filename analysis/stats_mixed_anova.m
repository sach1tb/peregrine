function stats_mixed_anova(type, cselect)


[bin_centers, units, jnk1, jnk2, id, id2]= get_options(type);

if nargin<2, cselect=1:numel(id); end

% anal='Variation in '; 
anal='Average ';



%%%%% STATS
data=csvread(['../data/data2_', type, '.csv']);
% reps=size(data,1)/numel(id);
reps=10;
data2=[];
for jj=cselect
    data2=[data2; data(reps*(jj-1)+1:reps*jj,:)];
end

% data2(isnan(data2))=rand;

X=data2(:);

% conditions
bsf=kron(ones(5,1), kron((1:numel(cselect))', ones(reps,1)));

% time
wsf=kron((1:5)',ones(numel(cselect)*reps,1));

% subject codes
sc=kron(ones(5,1),(1:numel(cselect)*reps)');

[SSQs, DFs, MSQs, Fs, Ps]=mixed_between_within_anova([X bsf, wsf, sc]);

fprintf('Mixed ANOVA, use only for *time-effect*\n\n');

% fprintf('Rows are conditions columns is time\n');
% 
% % p=anova2(data2,reps);
% pause()
Ps=cell2mat(Ps);
fprintf('Time effect (within subjects)=%.5f\nCondition effect(between subjects)=%.5f\nInteraction=%.5f\n', ...
                Ps(2), Ps(1), Ps(3));            

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