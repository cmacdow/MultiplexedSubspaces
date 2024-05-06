function [ndim,d,rPEV] = LocalDimRRR(src,trg,tidx,thresh)
%Camden M. - timeless
%estimates local dimensionality 
%splits local neural population into two spatially seperate
%populations and xval predicts activity in each and returns the #dim

%ToDo: 
%camden make sure not doing too mang remeaning of stuff. 
%confirm that they are organized by depth

if nargin <3; tidx = []; end
if nargin<4; thresh=0.8; end

if isempty(tidx) %just vanilla RRR
    [cvl_rrr,cvl_ridge, ~, ~] = briefRRR(src,trg,0); %this is just the local function below. Follows all the other RRR. 
    d = nanmean(1-cvl_rrr)/nanmean(bestLambda(cvl_ridge));
    ndim = find(d>=thresh,1,'first');
else %run on split trials to compute 'reliable variance    
    pev = RRR_reliable(src,trg,tidx,1);    %set last input to 0 to keep psth
   
    d = pev; %save off even the non-predicative dimensions
    
    %Get the predictive dimensions
    %at higher dimensions, it starts fitting data very poorly (e.g., negative
    %PEV, because it just noise. 
    pev(find(pev<0,1,'first'):end)=[];
    
    %total reliable variance
    rPEV = sum(pev);

    %get the num of dimensions      
    pev = cumsum(pev)/max(cumsum(pev));
    ndim = find(pev>=thresh,1,'first');

end


end

function [cvl_rrr, cvl_ridge,x, y] = briefRRR(x,y,rmvPSTH)
if nargin <3; rmvPSTH=1; end %remove the trial average activity

if rmvPSTH ==0 %don't subtract the psth
    %concatentate across trials 
    x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
    y = reshape(y,[size(y,1),size(y,2)*size(y,3)])'; 

    %mean center per neuron
    x = x-nanmean(x,1);
    y = y-nanmean(y,1);   

else %subtract psth (default)
    x = x-nanmean(x,3);
    y = y-nanmean(y,3);
    
    %concatentate across trials 
    x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
    y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';
end

% Number of cross validation folds.
cvNumFolds = 10;
lossMeasure = 'NSE'; % NSE stands for Normalized Squared Error
cvOptions = statset('crossval'); % Initialize default options for cross-validation.
numDimsUsedForPrediction = 1:min(30,size(x,2)); %prevent from going over num dim with small num neurons

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


function pev = RRR_reliable(x,y,tidx,rmvPSTH,ndim)
if nargin <4; rmvPSTH=1; end %remove the trial average activity
if nargin <5; ndim = min(size(y,1),30); end %keep consistent with rest of paper

if rmvPSTH ==0 %don't subtract the psth
    %concatentate across trials 
    x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
    y = reshape(y,[size(y,1),size(y,2)*size(y,3)])'; 

    %mean center per neuron
    x = x-nanmean(x,1);
    y = y-nanmean(y,1);   

else %subtract psth (default)
    x = x-nanmean(x,3);
    y = y-nanmean(y,3);
    
    %concatentate across trials 
    x = reshape(x,[size(x,1),size(x,2)*size(x,3)])';
    y = reshape(y,[size(y,1),size(y,2)*size(y,3)])';
end



%fit the model
[~,rrr_B,rrr_V] = ReducedRankRegress(y(tidx==0,:), x(tidx==0,:), ndim,'scale',false,'RIDGEINIT',0);

%Use fit model to predict y activity during withheld timepoints
ytest = y(tidx==1,:);
xtest = x(tidx==1,:);
pev = NaN(1,ndim);
% close all; figure;
for i = 1:ndim
    yhat = xtest * rrr_B(:,i) * rrr_V(:,i)';  
%     match intercepts
    yhat = yhat+(nanmean(ytest,1)-nanmean(yhat,1));
    pev(i) = (PercentExplainedVariance(ytest,yhat,0,0));
%     close all; figure; hold on; subplot(121); imagesc(ytest); colorbar; subplot(122); imagesc(yhat); colorbar; title(sprintf('dim %d pev%0.2g',i,pev(i))); pause();
end


%CAMDEN if intercepts are matched, then this must be true: 2∗∑(y−y^)(y^−y¯)=0
% find(2*sum((ytest-yhat).*(yhat-nanmean(ytest,1)))>1e-5);

end




% %     subplot(211); plot(ytest(:,1),yhat(:,1),'.'); title(sprintf('dim%d',i)); pause();
%     B = getB(ytest, xtest,rrr_B,rrr_V,i);
%     [~, yhat] = RegressPredict(ytest, xtest, B);
%     temp = (PercentExplainedVariance(ytest,yhat,0,0));
%     if i ==1
%         pev(2,i)=temp;
%     else
%         pev(2,i) = temp-sum(pev(2,1:i-1));
%     end
% %     subplot(212); plot(ytest(:,1),yhat(:,1),'.'); title(sprintf('dim%d',i)); pause();
% 
% function B = getB(Y,X,B_,V,dim)
% 
% [~, K] = size(Y);
% p = size(X, 2);
% 
% m = mean(X,1);
% 
% Bfull = B_/V;
% 
% if dim(1) == 0
%     B = zeros(p, K);
% else
%     B = Bfull*V( :, 1:dim(1) )*V( :, 1:dim(1) )';
% end
% 
% B = [ repmat( mean(Y,1), [1 numel(dim)] )-m*B; B ];
% 
% end

%%if you want to add ridge regression
% % % Number of cross validation folds.
% cvNumFolds = 5;
% lossMeasure = 'NSE'; % NSE stands for Normalized Squared Error
% cvOptions = statset('crossval'); % Initialize default options for cross-validation.
% 
% % %get the cross-validated value of lambda
% dMaxShrink = .5:.01:1;
% lambda = GetRidgeLambda(dMaxShrink, x,'scale',false);
% 
% cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
% 	(@RidgeRegress, Ytrain, Xtrain, Ytest, Xtest, lambda, ...
% 	'LossMeasure', lossMeasure,'scale',false); 
% 
% % Cross-validation routine.
% cvl_ridge = crossval(cvFun, y, x, ...
% 	  'KFold', cvNumFolds, ...
% 	'Options', cvOptions);
% 
% [~,idx] = bestLambda(cvl_ridge);

%Train model to predict Y on half the trials. 
% [~,rrr_B,rrr_V] = ReducedRankRegress(y(tidx==0,:), x(tidx==0,:), ndim,'scale',false,'RIDGEINIT',lambda(idx));










