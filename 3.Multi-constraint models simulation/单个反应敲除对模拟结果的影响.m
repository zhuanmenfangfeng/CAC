%% 先统计所有的反应条目：
rxn_names = {}
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建\Multiscale_Model_ensemble
for i = 1:500
    load(['model-',num2str(i)]);
    model = Model_thermo_;
    rxn_names = [rxn_names;model.rxns];
    clc
    sprintf('%d',i)
end
rxn_names = unique(rxn_names); 

% 删除生物质反应
rxn_names_1 = rxn_names;
r = find(contains(rxn_names_1,'biomass'));
rxn_names_1(r,:) = [];                            % 这些反应从待选择选择的集合中删除掉

% 合并可逆反应
f = @(x)x(3:end);
rxn_names_1 = cellfun(f,rxn_names_1,'UniformOutput',false);
rxn_names_1 = unique(rxn_names_1)                     % 数量是1533，正确了

% 删除交换反应
r = find(contains(rxn_names_1,'EX_'))
rxn_names_1(r,:) = [];

% 建立比生长速率模拟的矩阵
growth_mat = []
his_production = []

% 建立单反应敲除多约束模型求解算法
for i = 1:500
    cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建\Multiscale_Model_ensemble
    load(['model-',num2str(i)]);
    model = Model_thermo_;
    cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\3.多约束模型集模拟
    for ii = 1:length(rxn_names_1)
        model1 = model;                                                    % 删除相应反应采用什么方法呢？
        r = find(contains(model1.varNames,rxn_names_1{ii,1}));
        model1.A(:,r) = [];                                                % 模型反应条目删除工作
        model1.varNames(r,:) = [];
        model1.vartypes(r,:) = [];
        model1.var_lb(r,:) = [];
        model1.var_ub(r,:) = [];
        model1.f(r,:) = [];
        solution = gurobi_solve_CY(model1,'max');                          % 求解
        growth_mat(i,ii) = solution.objval;
    end
    clc
    sprintf('正在处理模型%d',i)
end

%%%%%%%%%% 这次同时观察两个东西，一个是统计比生长速率，一个是组氨酸的合成速率：而且能不能把这个东西加个速啊，我还有好几张图想放在补充材料里面了

rxn_names = {}
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建\Multiscale_Model_ensemble
for i = 1:500
    load(['model-',num2str(i)]);
    model = Model_thermo_;
    rxn_names = [rxn_names;model.rxns];
    clc
    sprintf('%d',i)
end
rxn_names = unique(rxn_names); 
































%% KNN+RF
% 首先登录上模型看看500个解分别对应的是哪些东西？
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\3.多约束模型集模拟
load rxn_names.mat

cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建\Multiscale_Model_ensemble
rxn_names_rxns_pos = zeros(500,length(rxn_names_1));
for i = 1:500
    load(['model-',num2str(i)]);
    model = Model_thermo_;
    for ii = 1:length(model.rxns)
        c = find(strcmp(model.rxns{ii}(3:end),rxn_names_1));           % 得把前缀拿走，要不然搜索不到相应的内容  
        rxn_names_rxns_pos(i,c) = 1;
    end
    clc
    sprintf('%d',i)
end
rxn_names_rxns_pos_1 = unique(rxn_names_rxns_pos,'rows');              % 里面检查到有一个内容是重复的，检查一下是哪个，这样就不用从

for i = 1:500
    if sum(sum(rxn_names_rxns_pos(i,:)' == rxn_names_rxns_pos')==1018)>1
        sprintf('%d',i);
        break
    end
end

%% 折叠特征
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\3.多约束模型集模拟
load simulation_1.mat

% 折叠无用的特征                （先折叠全是1或者全是0的特征，这步骤我也可以就是）
growth_mat = growth_mat>0.05
xlswrite('折叠前的特征集.xlsx',growth_mat,'Sheet1')
growth_mat_1 = sum(growth_mat)
c = find((500>growth_mat_1)&(growth_mat_1>0))
growth_mat_2 = growth_mat(:,c)

% 顺便对反应ID进行处理：
load rxn_names.mat
rxn_names_1 = rxn_names_1'
rxn_names_1 = rxn_names_1(:,c)
checkrxns = rxn_names_1
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
% 折叠完只有7个，但是也是可以用的,用unique函数能看出来确实是7个

% 保存折叠的特征：
growth_mat_2 = growth_mat_2'
growth_mat_3 = []
for i = 1:size(Res_rxns,1)
    c = find(strcmp(Res_rxns{i,1},rxn_names_1));
    growth_mat_3 = [growth_mat_3,growth_mat_2(:,c)];
end
save collapsed_essentiality.mat Res_rxns growth_mat_3

%% 开始聚类和训练
% Calculate the gene knockout distance between ensemble members, perform PCOA, 
% assign clusters with k-means, then classify with random forest.
% 计算模型集之间的组基因敲除距离，进行PCOA分析。用K均值进行聚类。

load collapsed_essentiality.mat

% 先生成一个对应的表格
Res = {}
for i = 1:size(growth_mat_3,1)
    for ii = 1:size(growth_mat_3,2)
        if growth_mat_3(i,ii)==1
            Res{i,ii} = 'True';
        else
            Res{i,ii} = 'False';
        end
    end
end
Res_rxns_names = Res_rxns(:,1);
Res_rxns_names = Res_rxns_names'
Res = [Res_rxns_names;Res]
 
Model_name_column = {}
for i = 2:size(Res,1)
    Model_name_column{i,1} = ['Model_',num2str(i-1)];
end
Res = [Model_name_column,Res]
xlswrite('模拟结果——折叠.xlsx',Res,'Sheet1')

%% 聚类完直接训练就行，不用考虑那么复杂
load collapsed_essentiality.mat

df_for_pcoa = growth_mat_3

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

% 模型训练过程以及参数选择可以用k-cv的方法以及模型泛化能力的测试，反正感觉这部分可以先跳过去了。
% 还有一个聚类系数，搞完就行了。
% 聚类系数：

Res_1 = []
Res_2 = []
for i = 1:length(idx)
    if idx(i,1)==1
        Res_1 = [Res_1;df_for_pcoa(i,:)];
    else
        Res_2 = [Res_2;df_for_pcoa(i,:)];
    end
end

Res_cluster = []
for i = 1:7
    Res_cluster(i,1) = sum(Res_1(:,i))./size(Res_1,1);
    Res_cluster(i,2) = sum(Res_2(:,i))./size(Res_2,1);
end

% 从上到下求解聚类系数：
for i = 1:size(Res_cluster,1)
    Res_cluster(i,3) = 1-min(Res_cluster(i,1:2))./max(Res_cluster(i,1:2));
end

% 随机森林参数优化：
df_for_pcoa = growth_mat_3

% 进行kmeans聚类，
opts = statset('Display','final');
[idx,C] = kmeans(df_for_pcoa,2,'Distance','cityblock','Replicates',5,'Options',opts)

% 将df_for_pcoa保存为pdFrame文件
xlswrite('折叠后的特征集',growth_mat_3,'Sheet1')

% 先生成一个对应的表格
Res = {}
for i = 1:size(growth_mat_3,1)
    for ii = 1:size(growth_mat_3,2)
        if growth_mat_3(i,ii)==1
            Res{i,ii} = 'True';
        else
            Res{i,ii} = 'False';
        end
    end
end
Res_rxns_names = Res_rxns(:,1);
Res_rxns_names = Res_rxns_names'
Res = [Res_rxns_names;Res]
 
Model_name_column = {}
for i = 2:size(Res,1)
    Model_name_column{i,1} = ['Model_',num2str(i-1)];
end

Res = [Model_name_column,Res]

%% 后续发现进行PcoA的意义其实没那么大，但是也可以先把PcoA分析给做出来：
py.importlib.import_module('skbio.stats.ordination')
py.importlib.import_module('skbio.stats.distance')

% 构建PdFrame结构
py.importlib.import_module('pandas');
pcoa_input = DissimilarityMatrix(squareform(pdist(df_for_pcoa,metric="hamming")),ids=df_for_pcoa.index)
pcoa_result = pcoa(pcoa_input, number_of_dimensions=3)

xlswrite('模拟结果.xlsx',growth_mat,'Sheet1')

[idx,C] = kmeans(df_for_pcoa,2,'Distance','cityblock','Replicates',5,'Options',opts)
xlswrite('聚类结果.xlsx',idx,'Sheet1')

%% 感觉还是用那个袋外交叉验证的方法。

c = cvpartition(10,'KFold',10)
TrainIndex = training(c,4)
TestIndex = training(c,2)

%% 贝叶斯优化算法的构建：
%%
maxMinLS = 20;
numTREE = optimizableVariable('numTREE',[1,500],'Type','integer');
minLS = optimizableVariable('minLS',[1,maxMinLS],'Type','integer');

hyperparametersRF = [numTREE;minLS];
results = bayesopt(@(params)oobErrRF(params,X,idx),hyperparametersRF,...
    'AcquisitionFunctionName','expected-improvement-plus','Verbose',0);

%% 





