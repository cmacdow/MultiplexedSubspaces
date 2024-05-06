function [cvl_ridge,rrr_B,rrr_V,cvl_rrr] = RRR_simple(x,y,lambda,xvalflag)
if nargin <3; lambda = []; end %refit ridge lambda
if nargin <3; xvalflag = 0; end %use cross validation (if you need rsq)
if size(x,3)>1
    %subtract the psth
    x = x-nanmean(x,3);
    y = y-nanmean(y,3);

    %concatentate across trials and pca
    x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
    y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';
end
% Number of cross validation folds.
cvNumFolds = 10;
lossMeasure = 'NSE'; % NSE stands for Normalized Squared Error
cvOptions = statset('crossval'); % Initialize default options for cross-validation.
numDimsUsedForPrediction = 1:30;

%% Full model | Ridge regression
if isempty(lambda)
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
    lambda = lambda(idx);
else
    cvl_ridge = [];
end

if xvalflag == 1
    rng('default')
    d = 1:min(numDimsUsedForPrediction(end),size(y,2));
    cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
        (@ReducedRankRegress, Ytrain, Xtrain, Ytest, Xtest, ...
        d, 'LossMeasure', lossMeasure,'scale',false,'RIDGEINIT',lambda);

    % Cross-validation routine.
    cvl_rrr = crossval(cvFun, y, x, ...
          'KFold', cvNumFolds, ...
        'Options', cvOptions);

    rrr_B=[];    
    rrr_V = [];
    
else
    %% Reduced Rank Regression 
    rng('default')
    d = 1:min(numDimsUsedForPrediction(end),size(y,2));
    % Reduced Rank Regression Full with the optimal dimensionality
    [~,rrr_B,rrr_V] = ReducedRankRegress(y, x, d,'scale',false,'RIDGEINIT',lambda);
    cvl_rrr =[];
end

end %function end












