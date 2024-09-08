% 我不能完全折叠成一样的特征
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\3.多约束模型集模拟
load simulation_1.mat
load rxn_names.mat

% 先折叠完全相同的特征：
growth_mat_2 = []
rxn_names_2 = {}
for i = 1:size(growth_mat,2)
    if length(unique(growth_mat(:,i)))>1
       growth_mat_2 = [growth_mat_2,growth_mat(:,i)];               % 确实少了一些
       rxn_names_2 = [rxn_names_2,rxn_names_1(i,1)];
    end
end

% 前面是特征的折叠，这个其实还是折叠特征了，只不过现在是特征间的相似度，上面是根据特征内的相似度。

checkrxns = rxn_names_2
usedrxns = {}

sprintf("Genes before collapse: " + num2str(length(checkrxns)))
Res_rxns = {}
b = 1

for i = 1:length(checkrxns)
    rxn1 = checkrxns(1,i);
    if sum(strcmp(rxn1,usedrxns))==0
        identical_rxns = rxn1;
        usedrxns = [usedrxns;rxn1];
        for ii = 1:length(checkrxns)
            rxn2 = checkrxns{1,ii}; 
            if ~strcmp(rxn1,rxn2)&sum(strcmp(rxn2,usedrxns))==0
                if sum(growth_mat_2(:,i)==growth_mat_2(:,ii))==500;
                    identical_rxns = [identical_rxns;rxn2];
                    usedrxns = [usedrxns;rxn2];
                end
            end
        end
        for ii = 1:length(identical_rxns)
            Res_rxns{b,ii} = identical_rxns{ii,1};
        end
        b = b+1;
    end
end
sprintf("Genes after collapse: " + num2str(size(Res_rxns,2)));

% 我擦 折叠完还有400个特征:

% 保存折叠的特征：
growth_mat_2 = growth_mat_2'
growth_mat_3 = [];

for i = 1:size(Res_rxns,1)
    c = find(strcmp(Res_rxns{i,1},rxn_names_2));
    growth_mat_3 = [growth_mat_3,growth_mat_2(:,c)];
end

cd
save collapsed_essentiality_for_GBDT.mat Res_rxns growth_mat_3
%% 直接进行PCA分析，不用考虑的太复杂

load collapsed_essentiality_for_GBDT.mat
load ensmeble_rxn.mat

% 统计每个模型集中反应出现频率 特征部分：

total_rxn = {}
for i = 1:length(Res_)
    total_rxn = [total_rxn;Res_{i,1}];
end
total_rxn = unique(total_rxn)
total_rxn = total_rxn'

train_data = zeros(500,length(total_rxn))
for i = 1:length(Res_)
    data = Res_{i,1};
    for ii = 1:length(data)
        [~,c] = find(strcmp(data{ii,1},total_rxn));
        train_data(i,c) = 1;
    end
end



% 我得先进行PCA分析：
[coeff,score,latent] = pca(growth_mat_3);                    % 这部分的图可以放在补充材料中。
idx = score(:,1);

% 随机森林建立模型：
trees = 500;                                                    % 决策树数目
leaf  = 5;                                                      % 最小叶子数
OOBPrediction = 'on';                                           % 打开误差图
OOBPredictorImportance = 'on';                                  % 计算特征重要性
Method = 'regression';                                          % 分类还是回归
net = TreeBagger(trees, train_data, idx, 'OOBPredictorImportance', OOBPredictorImportance,...
      'Method', Method, 'OOBPrediction', OOBPrediction, 'minleaf', leaf);
importance = net.OOBPermutedPredictorDeltaError;  % 重要性

plot(importance) % 这是这种方法获得的重要性顺序 还有 之前的方法的重要性排序。 我可以加一张参数优化图。

xlswrite('重要性结果.xlsx',importance,'Sheet1')

clear

% 再计算一个cluster_ratio
cluster_ratio = []

total_rxn = total_rxn'
for i = 1:length(total_rxn)
    r_contains = [];
    r_nocontains = [];
    for ii = 1:length(Res_)
        if sum(strcmp(total_rxn{i,1},Res_{ii,1}))>0
            r_contains = [r_contains;ii];
        else
            r_nocontains = [r_nocontains;ii];
        end
    end
    fr1 = mean(idx(r_contains,:));
    fr2 = mean(idx(r_nocontains,:));
    total_rxn{i,2} = 1-min([fr1,fr2])./max([fr1,fr2]);
end

xlswrite('重要性结果.xlsx',total_rxn,'new_method_cluste_ratio')

%% 传统方法对比生长速率重要性的
% 传统方法，对结果进行聚类。
clear

cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\3.多约束模型集模拟
load collapsed_essentiality.mat
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\8.模型集模拟策展方法
load ensmeble_rxn.mat

total_rxn = {}
for i = 1:length(Res_)
    total_rxn = [total_rxn;Res_{i,1}];
end
total_rxn = unique(total_rxn)
total_rxn = total_rxn'

train_data = zeros(500,length(total_rxn))
for i = 1:length(Res_)
    data = Res_{i,1};
    for ii = 1:length(data)
        [~,c] = find(strcmp(data{ii,1},total_rxn));
        train_data(i,c) = 1;
    end
end
df_for_pcoa = train_data

% 进行kmeans聚类，
opts = statset('Display','final');
[idx,C] = kmeans(df_for_pcoa,2,'Distance','cityblock','Replicates',5,'Options',opts)

% 随机森林建立模型：
trees = 500;                                      % 决策树数目
leaf  = 5;                                        % 最小叶子数
OOBPrediction = 'on';                             % 打开误差图
OOBPredictorImportance = 'on';                    % 计算特征重要性
Method = 'regression';                            % 分类还是回归
net = TreeBagger(trees, df_for_pcoa, idx, 'OOBPredictorImportance', OOBPredictorImportance,...
      'Method', Method, 'OOBPrediction', OOBPrediction, 'minleaf', leaf);
importance = net.OOBPermutedPredictorDeltaError;  % 重要性
xlswrite('重要性结果.xlsx',importance,'传统方法二元变量')
xlswrite('重要性结果.xlsx',total_rxn(:,end-83:end)','reaction_importance')

% 
r1 = find(idx==1)
r2 = find(idx==2)
total_rxn = total_rxn'

data_1 = Res_(r1,:)
data_2 = Res_(r2,:)


for i = 1:length(total_rxn)
    b = 0;
    for ii = 1:length(data_1)
        if sum(strcmp(total_rxn{i,1},data_1{ii,1}))>0
           b = b+1;
        end
    end
    fr1 = b./length(data_1);

    b = 0;
    for ii = 1:length(data_2)
        if sum(strcmp(total_rxn{i,1},data_2{ii,1}))>0
            b = b+1;
        end
    end
    fr2 = b./length(data_1);
    total_rxn{i,2} = 1-min([fr1,fr2])./max([fr1,fr2]);

end
xlswrite('重要性结果.xlsx',total_rxn,'old_method_cluste_ratio');

% 归一化
[Y,PS] = mapminmax(X',0,1)
Y = Y'






