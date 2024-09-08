load collapsed_essentiality.mat
X = growth_mat_3
maxMinLS = 20;
numTREE = optimizableVariable('numTREE',[1,500],'Type','integer');
minLS = optimizableVariable('minLS',[1,maxMinLS],'Type','integer');

% 进行kmeans聚类，
opts = statset('Display','final');
[idx,C] = kmeans(X,2,'Distance','cityblock','Replicates',5,'Options',opts)
hyperparametersRF = [numTREE;minLS];
results = bayesopt(@(params)oobErrRF(params,X,idx),hyperparametersRF,'AcquisitionFunctionName','expected-improvement-plus','Verbose',0);

function oobErr = oobErrRF(params,X,idx)
%oobErrRF Trains random forest and estimates out-of-bag quantile error
%   oobErr trains a random forest of 300 regression trees using the
%   predictor data in X and the parameter specification in params, and then
%   returns the out-of-bag quantile error based on the median. X is a table
%   and params is an array of OptimizableVariable objects corresponding to
%   the minimum leaf size and number of predictors to sample at each node.
randomForest = TreeBagger(params.numTREE,X,idx,'Method','regression',...
    'OOBPrediction','on','MinLeafSize',params.minLS,...
    'OOBPrediction', 'on');
oobErr = oobQuantileError(randomForest);
end

bestOOBErr = results.MinObjective                        % 最好的目标函数值
bestHyperparameters = results.XAtMinObjective






