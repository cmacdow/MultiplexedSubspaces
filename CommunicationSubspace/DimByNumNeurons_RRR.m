function DimByNumNeurons_RRR(data,normtype)

%Camden - timeless
%combines data across recordings and motifs to return a single variable

if nargin <2; normtype = 'mean'; end
all_areas = data{1}(1).area_label;


idx = strcmp(data{cur_rec}(cur_motif).area_label,all_areas{cur_area});
fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
fr = normalizeToBaseline(fr,1:2,normtype);
src = fr(:,3:end,:); %source
fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas
fr = normalizeToBaseline(fr,1:2,normtype);
trg = fr(:,3:end,:); %source


%% run the loop
%split number of targets into divisions of ten
n = floor(linspace(1,size(trg,1),11));
n = n(2:end);
rng('default')
d = NaN(numel(n),30);
for i = 1:numel(n)  
    fprintf('\nworking on %d or %d',i,numel(n))
   [cvl_rrr, cvl_ridge] = briefRRR(src,trg(randperm(n(end),n(i)),:,:));
   d(i,:) = nanmean(1-cvl_rrr)/nanmean(bestLambda(cvl_ridge,0));
end


%% varying source
%split number of targets into divisions of ten
n = floor(linspace(1,size(src,1),11));
n = n(2:end);
rng('default')
d = NaN(numel(n),30);
for i = 1:numel(n)  
    fprintf('\nworking on %d or %d',i,numel(n))
   [cvl_rrr, cvl_ridge] = briefRRR(src(randperm(n(end),n(i)),:,:),trg);
   d(i,:) = nanmean(1-cvl_rrr)/nanmean(bestLambda(cvl_ridge,0));
end




% 
% x = NaN(numel(data),size(data{1},2),numel(all_areas),2000);
% for cur_rec = 1:numel(data) %loop through each recording
%     for cur_motif = 1:size(data{cur_rec},2)%loop throuch each motif
%         for cur_area = 1:numel(all_areas)            
%             idx = strcmp(data{cur_rec}(cur_motif).area_label,all_areas{cur_area});
%             if sum(idx) >0               
%                 fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); %all other areas, i.e., to match the betas
%                 fr = normalizeToBaseline(fr,1:2,normtype);
%                 src = fr(:,3:end,:); %source
%                 fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); %all other areas, i.e., to match the betas
%                 fr = normalizeToBaseline(fr,1:2,normtype);
%                 trg = fr(:,3:end,:); %source
%                 [cvl_rrr, cvl_ridge] = briefRRR(src,trg);
%             end
%         end
%     end
% end


end %function end


function [cvl_rrr, cvl_ridge] = briefRRR(x,y)

%subtract the psth
x = x-nanmean(x,3);
y = y-nanmean(y,3);

%concatentate across trials 
x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';

% Number of cross validation folds.
cvNumFolds = 10;
lossMeasure = 'NSE'; % NSE stands for Normalized Squared Error
cvOptions = statset('crossval'); % Initialize default options for cross-validation.
numDimsUsedForPrediction = 1:30;

%% Full model | Ridge regression
rng('default')
dMaxShrink = .5:.01:1;
lambda = GetRidgeLambda(dMaxShrink, x,'scale',false);

cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
	(@RidgeRegress, Ytrain, Xtrain, Ytest, Xtest, lambda, ...
	'LossMeasure', lossMeasure,'scale',false); 

% Cross-validation routine.
cvl_ridge = crossval(cvFun, y, x, ...
	  'KFold', cvNumFolds, ...
	'Options', cvOptions);
[~,idx] = bestLambda(cvl_ridge);
%% Reduced Rank Regression XVal
rng('default')
d = 1:min(numDimsUsedForPrediction(end),size(y,2));
cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
	(@ReducedRankRegress, Ytrain, Xtrain, Ytest, Xtest, ...
	d, 'LossMeasure', lossMeasure,'scale',false,'RIDGEINIT',lambda(idx));

% Cross-validation routine.
cvl_rrr = crossval(cvFun, y, x, ...
	  'KFold', cvNumFolds, ...
	'Options', cvOptions);
end

