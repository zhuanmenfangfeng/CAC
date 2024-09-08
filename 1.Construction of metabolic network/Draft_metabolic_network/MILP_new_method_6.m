load Unversal-model_icw773metadd_rxnadd_rxnrep.mat;
Model = Universal_Model;
[~,~,block_reaction] = xlsread('全局模型中的block_reaction.xlsx','Sheet1');
[~,~,repeat_reaction] = xlsread('全局反应中重复的反应.xlsx','Sheet1');
Model = removeRxns(Model,block_reaction);
Model = removeRxns(Model,repeat_reaction);                                 % ok，现在变成了8000多个反应，

% 限制全部的交换反应：
r_EX_ = find(contains(Model.rxns,'EX_'));
Model.lb(r_EX_,:) = 0;

% 增加葡萄糖耗散反应：
Model = addReaction(Model,'Glucose_dissipation_reaction','metaboliteList',{'cpd00027[c0]'},'stoichCoeffList',[-1], 'reversible',false);
Model = changeObjective(Model,'Glucose_dissipation_reaction');

% 求解
solution = optimizeCbModel(Model)
model2xls(Model,'模型240118.xlsx')
%% 处理表格
[~,~,data] = xlsread('data中仅保留1329号反应.xlsx','Reaction List')
r = []
for i = 1:length(data)
    if contains(data{i,3},'cpd00027[c0]') & data{i,9} ~= 0
       r = [r;i];
    end
end
data1 = data(r,:)
xlswrite('模型240118.xlsx',data1,'通量流向葡萄糖的方程')

%% 直接在此基础上移除一下之前确定的反应条目：（rxn05746_c0、rxn12548_c0、rxn90123_c0）
clear;
load Unversal-model_icw773metadd_rxnadd_rxnrep.mat;
Model = Universal_Model;
[~,~,block_reaction] = xlsread('全局模型中的block_reaction.xlsx','Sheet1');
[~,~,repeat_reaction] = xlsread('全局反应中重复的反应.xlsx','Sheet1');
Model = removeRxns(Model,block_reaction);
Model = removeRxns(Model,repeat_reaction);                                 % ok，现在变成了8000多个反应，
wait_remove = {'rxn05746_c0';'rxn12548_c0';'rxn90123_c0'}
Model = removeRxns(Model,wait_remove)

% 限制营养条件
r_EX_ = find(contains(Model.rxns,'EX_'));
Model.lb(r_EX_,:) = 0;

[~,~,basic_substance] = xlsread('biolog_base_composition.csv')
basic_substance = basic_substance(:,2)
basic_substance(1,:) = []

% 由于icw773生物量反应的特殊性，需要额外添加这些物质
waitaddmets_new = {'cpd00244[e0]';'cpd00104[e0]';'cpd00099[e0]';'cpd00027[e0]'}
wait_open_rxn = [waitaddmets_new;basic_substance]
f = @(x)strrep(x,'_e','[e0]')
wait_open_rxn = cellfun(f,wait_open_rxn,'UniformOutput',false)

% 确定对应的交换反应
for i = 1:length(wait_open_rxn)
    r = find(strcmp(wait_open_rxn{i,1},Universal_Model.mets));
    c = find(Universal_Model.S(r,:)==-1);
    for ii = 1:length(c)
        if length(find(Universal_Model.S(:,c(ii))~=0))==1;
            wait_open_rxn{i,2} = Universal_Model.rxns{c(ii),1};
        end
    end
end

wait_open_rxn(4,:) = []

% 设置相应的反应速度
for i = 1:length(wait_open_rxn)
    c_ = find(strcmp(Model.rxns,wait_open_rxn{i,2}));
    Model.lb(c_,1) = -100;
end

Model = changeObjective(Model,'add_CG_biomass cgl ATCC13032')
solution = optimizeCbModel(Model)                          % 可以 确实没有比生长速率

% 打开葡萄糖交换反应
[~,~,basic_substance] = xlsread('biolog_base_composition.csv')
basic_substance = basic_substance(:,2)
basic_substance(1,:) = []

% 由于icw773生物量反应的特殊性，需要额外添加这些物质
waitaddmets_new = {'cpd00244[e0]';'cpd00104[e0]';'cpd00099[e0]';'cpd00027[e0]'}
wait_open_rxn = [waitaddmets_new;basic_substance]
f = @(x)strrep(x,'_e','[e0]')
wait_open_rxn = cellfun(f,wait_open_rxn,'UniformOutput',false)

% 确定对应的交换反应
for i = 1:length(wait_open_rxn)
    r = find(strcmp(wait_open_rxn{i,1},Universal_Model.mets));
    c = find(Universal_Model.S(r,:)==-1);
    for ii = 1:length(c)
        if length(find(Universal_Model.S(:,c(ii))~=0))==1;
            wait_open_rxn{i,2} = Universal_Model.rxns{c(ii),1};
        end
    end
end
                                                           
% 设置相应的反应速度
for i = 1:length(wait_open_rxn)
    c_ = find(strcmp(Model.rxns,wait_open_rxn{i,2}));
    Model.lb(c_,1) = -100;
end

Model = changeObjective(Model,'add_CG_biomass cgl ATCC13032')
solution = optimizeCbModel(Model)                                          % 可以 这时候又有比生长速率了，那就可以进行正常的gap_filling模块。

%% 再检查一下他的碳源利用情况：
% 看其在其他碳源数据上的表现：

clear;
load Unversal-model_icw773metadd_rxnadd_rxnrep.mat;
Model = Universal_Model;
[~,~,block_reaction] = xlsread('全局模型中的block_reaction.xlsx','Sheet1');
[~,~,repeat_reaction] = xlsread('全局反应中重复的反应.xlsx','Sheet1');
Model = removeRxns(Model,block_reaction);
Model = removeRxns(Model,repeat_reaction);                                 % ok，现在变成了8000多个反应，
wait_remove = {'rxn05746_c0';'rxn12548_c0';'rxn90123_c0'}
Model = removeRxns(Model,wait_remove)

% 限制营养条件
r_EX_ = find(contains(Model.rxns,'EX_'));
Model.lb(r_EX_,:) = 0;

[~,~,basic_substance] = xlsread('biolog_base_composition.csv')
basic_substance = basic_substance(:,2)
basic_substance(1,:) = []

% 由于icw773生物量反应的特殊性，需要额外添加这些物质
waitaddmets_new = {'cpd00244[e0]';'cpd00104[e0]';'cpd00099[e0]';'cpd00027[e0]'}
wait_open_rxn = [waitaddmets_new;basic_substance]
f = @(x)strrep(x,'_e','[e0]')
wait_open_rxn = cellfun(f,wait_open_rxn,'UniformOutput',false)

% 确定对应的交换反应
for i = 1:length(wait_open_rxn)
    r = find(strcmp(wait_open_rxn{i,1},Universal_Model.mets));
    c = find(Universal_Model.S(r,:)==-1);
    for ii = 1:length(c)
        if length(find(Universal_Model.S(:,c(ii))~=0))==1;
            wait_open_rxn{i,2} = Universal_Model.rxns{c(ii),1};
        end
    end
end

wait_open_rxn(4,:) = []

% 设置相应的反应速度
for i = 1:length(wait_open_rxn)
    c_ = find(strcmp(Model.rxns,wait_open_rxn{i,2}));
    Model.lb(c_,1) = -100;
end

Model = changeObjective(Model,'add_CG_biomass cgl ATCC13032')
solution = optimizeCbModel(Model)                          % 可以 确实没有比生

% 查看在其他碳源上的利用能力：
[~,~,carbon_source_test] = xlsread('plata_thresholded.csv')
carbon_source_substance = split(carbon_source_test{1,1},'	')
carbon_source_substance(1,:) = []

f = @(x)strrep(x,'_e','[e0]')
carbon_source_substance = cellfun(f,carbon_source_substance,'UniformOutput',false)
for i = 1:length(carbon_source_substance)
    if sum(strcmp(carbon_source_substance{i,1},Universal_Model.mets))>0
        carbon_source_substance{i,2} = 0;
    else
        carbon_source_substance{i,2} = 1;
    end
end     

Res = {}
for i = 1:length(carbon_source_test)
    a = split(carbon_source_test{i,1},'	');
    a = a';
    Res = cat(1,Res,a);
end
r = find(strcmp(Res(:,1),'Corynebacterium glutamicum'))

Res(1,:) = cellfun(f,Res(1,:),'UniformOutput',false)
for i = 1:length(carbon_source_substance)
    c = find(strcmp(carbon_source_substance{i,1},Res(1,:)));
    carbon_source_substance{i,3} = Res{r,c};
end

r = find(cell2mat(carbon_source_substance(:,2))==1)
carbon_source_substance(r,:) = []

for i = 1:length(carbon_source_substance)
    r = find(strcmp(carbon_source_substance{i,1},Model.mets));
    c = find(Model.S(r,:)==-1);
    for ii = 1:length(c)
        if length(find(Model.S(:,c(ii))~=0))==1;
            carbon_source_substance{i,4} = Model.rxns{c(ii),1};
        end
    end
end

for i = 1:length(carbon_source_substance)
    Model_1 = Model;
    Model_1 = changeRxnBounds(Model_1,carbon_source_substance{i,4},-100,'l');
    solution = optimizeCbModel(Model_1);
    carbon_source_substance{i,5} = solution.f;
end                                                         % 好像看起来还可以，不管有没有碳源都能力利用

% 补gap用到哪些数据：
r1 = find(strcmp(carbon_source_substance(:,3),'True'));
r2 = find(cell2mat(carbon_source_substance(:,5))>0);
could_applicate = intersect(r1,r2);
carbon_source_substance_1 = carbon_source_substance(could_applicate,:)

% 刨除一些初始模型可利用的数据：
load need_fill_gap.mat
for i = 1:length(carbon_source_substance)
    if sum(strcmp(carbon_source_substance{i,1},carbon_source_substance_1(:,1)))>0
        r = find(strcmp(carbon_source_substance{i,1},carbon_source_substance_1(:,1)));
        carbon_source_substance{i,5} = carbon_source_substance_1{r,5};
    end
end
save finall_gap_fill.mat carbon_source_substance






