data_rxn = readtable('GCA_000196335.1_ASM19633v1_genomic-all-Reactions.txt');
data_transporter = readtable('GCA_000196335.1_ASM19633v1_genomic-Transporter.txt');

%% 给每个反应赋以相应的分值，为每个反应计算相应的权重：
% 将所有的rxn的bitscore的得分改成权重
r = find(isnan(data_rxn.bitscore))
data_rxn.bitscore(r,:) = 0

% high_evi_rxn_BS = 200
% dummy_weight = 100
% min_bs_for_core = 50
% data_rxn.bitscore_normalization = (data_rxn.bitscore - high_evi_rxn_BS) *
% ((0.005 - dummy_weight)/(high_evi_rxn_BS - min_bs_for_core)) + 0.005
% 不得是把所有的bitscore放在一起归一化么？
a = mapminmax(data_rxn.bitscore',0,100)
data_rxn.bitscore_normalization = a'

b = mapminmax(data_transporter.bitscore',0,100)
data_transporter.bitscore_normalization = b'
save recorde_rxn_score.mat data_transporter data_rxn

%% 查看全局模型在最小培养基上的碳源利用能力：

load Unversal-model_icw773metadd_rxnadd_rxnrep.mat

% 先把所有的矩阵整合在一起：

[~,~,carbon_source_test] = xlsread('plata_thresholded.csv')
[~,~,basic_substance] = xlsread('biolog_base_composition.csv')
carbon_source_substance = split(carbon_source_test{1,1},'	')
carbon_source_substance(1,:) = []
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

Universal_Model = changeObjective(Universal_Model,'add_CG_biomass cgl ATCC13032')
optimizeCbModel(Universal_Model)      % 比生长速率正常

% 关闭模型中的交换反应
r3 = find(contains(Universal_Model.rxns,'_EX_'))
Universal_Model = changeRxnBounds(Universal_Model, Universal_Model.rxns(r3,:), 0, 'l')
optimizeCbModel(Universal_Model)      % 比生长速率为0

% 打开模型中的设定好交换反应
wait_open_rxn(4,:) = []
wait_open_rxn(22,:) = []
wait_open_rxn(2,:) = []
for i = 1:size(wait_open_rxn,1)
    r4 = find(strcmp(Universal_Model.rxns,wait_open_rxn{i,2}));
    Universal_Model = changeRxnBounds(Universal_Model, Universal_Model.rxns{r4,1}, -100, 'l');
end
optimizeCbModel(Universal_Model)                             % 确实可以生长，比生长速率为6

% 在全局模型中只开这些交换反应：
load Unversal-model_icw773metadd_rxnadd_rxnrep.mat
r3 = find(contains(Universal_Model.rxns,'EX_'))
Universal_Model = changeRxnBounds(Universal_Model, Universal_Model.rxns(r3,:), 0, 'l')
Universal_Model = changeObjective(Universal_Model,'add_CG_biomass cgl ATCC13032')
optimizeCbModel(Universal_Model)      % 比生长速率为0

for i = 1:size(wait_open_rxn,1)
    r4 = find(strcmp(Universal_Model.rxns,wait_open_rxn{i,2}));
    Universal_Model = changeRxnBounds(Universal_Model, Universal_Model.rxns{r4,1}, -100, 'l');
end
optimizeCbModel(Universal_Model)                         % 比生长速率变成了12

% 看其在其他碳源数据上的表现：
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
    r = find(strcmp(carbon_source_substance{i,1},Universal_Model.mets));
    c = find(Universal_Model.S(r,:)==-1);
    for ii = 1:length(c)
        if length(find(Universal_Model.S(:,c(ii))~=0))==1;
            carbon_source_substance{i,4} = Universal_Model.rxns{c(ii),1};
        end
    end
end

load Unversal-model_icw773metadd_rxnadd_rxnrep.mat
r3 = find(contains(Universal_Model.rxns,'EX_'))
Universal_Model = changeRxnBounds(Universal_Model, Universal_Model.rxns(r3,:), 0, 'l')
Universal_Model = changeObjective(Universal_Model,'add_CG_biomass cgl ATCC13032')
wait_open_rxn(4,:) = []     % 删除葡萄糖交换反应

Universal_Model = changeRxnBounds(Universal_Model, wait_open_rxn(:,2), -100, 'l');

for i = 1:length(carbon_source_substance)
    Universal_Model_test = Universal_Model;
    Universal_Model_test  = changeRxnBounds(Universal_Model_test, carbon_source_substance{i,4}, -100, 'l');
    sol = optimizeCbModel(Universal_Model_test);
    carbon_source_substance{i,5} = sol.f;
end



%% 形成混合整数线性规划
% 在单个培养条件下尝试一下：
% 读取模型：
load Unversal-model_icw773metadd_rxnadd_rxnrep.mat
load recorde_rxn_score.mat
% 先把所有的矩阵整合在一起：

[~,~,carbon_source_test] = xlsread('plata_thresholded.csv')
[~,~,basic_substance] = xlsread('biolog_base_composition.csv')
carbon_source_substance = split(carbon_source_test{1,1},'	')
carbon_source_substance(1,:) = []
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

% 只打开初始模型中的交换反应：
r1 = find(~contains(Universal_Model.rxns,'rep_'))
r2 = find(~contains(Universal_Model.rxns,'add_'))
r = intersect(r1,r2)
Universal_Model = changeRxnBounds(Universal_Model, Universal_Model.rxns(r,:), 0, 'b')
Universal_Model = changeObjective(Universal_Model,'add_CG_biomass cgl ATCC13032')
optimizeCbModel(Universal_Model)      % 比生长速率正常

% 关闭模型中的交换反应
r3 = find(contains(Universal_Model.rxns,'_EX_'))
Universal_Model = changeRxnBounds(Universal_Model, Universal_Model.rxns(r3,:), 0, 'l')
optimizeCbModel(Universal_Model)      % 比生长速率为0

% 打开模型中的设定好交换反应
for i = 1:size(wait_open_rxn,1)
    r4 = find(strcmp(Universal_Model.rxns,wait_open_rxn{i,2}));
    Universal_Model = changeRxnBounds(Universal_Model, Universal_Model.rxns{r4,1}, -100, 'l');
end
optimizeCbModel(Universal_Model)                             % 确实可以生长，比生长速率为6

% 在全局模型中只开这些交换反应：
load Unversal-model_icw773metadd_rxnadd_rxnrep.mat
r3 = find(contains(Universal_Model.rxns,'EX_'))
Universal_Model = changeRxnBounds(Universal_Model, Universal_Model.rxns(r3,:), 0, 'l')
Universal_Model = changeObjective(Universal_Model,'add_CG_biomass cgl ATCC13032')
optimizeCbModel(Universal_Model)      % 比生长速率为0

for i = 1:size(wait_open_rxn,1)
    r4 = find(strcmp(Universal_Model.rxns,wait_open_rxn{i,2}));
    Universal_Model = changeRxnBounds(Universal_Model, Universal_Model.rxns{r4,1}, -100, 'l');
end
optimizeCbModel(Universal_Model)                         % 比生长速率变成了12

% 看其在其他碳源数据上的表现：
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
    r = find(strcmp(carbon_source_substance{i,1},Universal_Model.mets));
    c = find(Universal_Model.S(r,:)==-1);
    for ii = 1:length(c)
        if length(find(Universal_Model.S(:,c(ii))~=0))==1;
            carbon_source_substance{i,4} = Universal_Model.rxns{c(ii),1};
        end
    end
end

load Unversal-model_icw773metadd_rxnadd_rxnrep.mat
r3 = find(contains(Universal_Model.rxns,'EX_'))
Universal_Model = changeRxnBounds(Universal_Model, Universal_Model.rxns(r3,:), 0, 'l')
Universal_Model = changeObjective(Universal_Model,'add_CG_biomass cgl ATCC13032')
wait_open_rxn(4,:) = []     % 删除葡萄糖交换反应

Universal_Model = changeRxnBounds(Universal_Model, wait_open_rxn(:,2), -100, 'l');

for i = 1:length(carbon_source_substance)
    Universal_Model_test = Universal_Model;
    Universal_Model_test  = changeRxnBounds(Universal_Model_test, carbon_source_substance{i,4}, -100, 'l');
    sol = optimizeCbModel(Universal_Model_test);
    carbon_source_substance{i,5} = sol.f;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 上述内容为Carbon_source_utilized 中的code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 看到有多少 Ture 的行

r = find(strcmp(carbon_source_substance(:,3),'False'))
carbon_source_substance_1 = carbon_source_substance
carbon_source_substance_1(r,:) = []
nutrition_conditions = size(carbon_source_substance_1,1)
%% MILP 矩阵初步构建
% 先尝试1培养基的条件,这次反应得分要用 zscore的方法进行归一化处理，并且拆分方程的可逆性。
load Unversal-model_icw773metadd_rxnadd_rxnrep.mat

nutrition_conditions = 1;
Model.A = sparse(nutrition_conditions.*size(Universal_Model.S,1),nutrition_conditions.*size(Universal_Model.S,2));
Model.constraintType = {};
Model.constraintNames = {};
Model.rhs = [];
b = 1;
% 增加1种培养条件的约束
for i = 1:nutrition_conditions
    for ii=1:size(Universal_Model.mets,1)
        Model.constraintType{b,1} = '=';
        Model.constraintNames{b,1} = strcat(num2str(i),'_',Universal_Model.mets{ii});  % 代谢物 质量平衡约束
        Model.rhs(b,1) = 0;
        Model.A(b,[((i-1).*length(Universal_Model.rxns)+1):i.*length(Universal_Model.rxns)]) = Universal_Model.S(ii,:);
        b = b + 1;
    end
end

% 增加1种培养条件的变量
Model.var_lb = [];
Model.var_ub = [];
Model.varNames = {}
Model.vartypes = []
for i = 1:nutrition_conditions
    for ii=1:size(Universal_Model.rxns,1)
        Model.var_lb = [Model.var_lb;Universal_Model.lb(ii)];
        Model.var_ub = [Model.var_ub;Universal_Model.ub(ii)];
        Model.varNames = [Model.varNames;strcat(num2str(i),'_',Universal_Model.rxns{ii})];
        Model.vartypes = [Model.vartypes;'C'];
    end
    clc
    sprintf('已经进行了%d',i)
end
Model.f = zeros(size(Model.A,2),1);
Model.vartypes = num2cell(Model.vartypes)

% 增加二元变量，用来限制反应数量
r1 = find(~contains(Universal_Model.rxns,'EX_'));
r2 = find(~contains(Universal_Model.rxns,'biomass'));
r = intersect(r1,r2)
rxn_wait_selected = Universal_Model.rxns(r,:);
for i = 1:length(r)
    Model = addNewVariable_CY(Model, strcat('select_or_not',rxn_wait_selected{i},'_f'),'B',[0,1]);
    Model = addNewVariable_CY(Model, strcat('select_or_not',rxn_wait_selected{i},'_r'),'B',[0,1]);
    clc
    sprintf('已经进行了%d',i)
end

tic
for i = 1:nutrition_conditions
    for ii = 1:length(r)
        Model.constraintNames = [Model.constraintNames;strcat('Constraint_lb_',num2str(i),'_',Universal_Model.rxns{r(ii),1})];
        Model.constraintNames = [Model.constraintNames;strcat('Constraint_ub_',num2str(i),'_',Universal_Model.rxns{r(ii),1})];
        Model.constraintNames = [Model.constraintNames;strcat('Constraint_ub_lb',num2str(i),'_',Universal_Model.rxns{r(ii),1})];

        Model.constraintType = [Model.constraintType;{'>'}];
        Model.constraintType = [Model.constraintType;{'<'}];
        Model.constraintType = [Model.constraintType;{'<'}];

        Model.rhs = [Model.rhs;0;0;1];
    end
    clc
    sprintf('%d',i)
end
toc

% 大矩阵中填写相应的系数
Model_ = Model

tic
for i = 1:nutrition_conditions
    for ii = 1:length(r)
        c_ = find(strcmp(Model_.varNames,strcat(num2str(i),'_',Universal_Model.rxns{r(ii),1})));
        pos_f = find(strcmp(Model_.varNames,strcat('select_or_not',rxn_wait_selected{i},'_f')));
        pos_r = find(strcmp(Model_.varNames,strcat('select_or_not',rxn_wait_selected{i},'_r')));

        a = sparse(3,size(Model_.A,2));
        CLHS.pos = [c_,pos_f,pos_r];
        CLHS.coeffs_1 = [1,-1e-6,1e6];
        CLHS.coeffs_2 = [1,-1e6,1e-6];
        
        a(1,CLHS.pos) = CLHS.coeffs_1;
        a(2,CLHS.pos) = CLHS.coeffs_2;
        a(3,[pos_f,pos_r]) = [1,1];
        
        Model_.A = [Model_.A;a];
    end
    clc
    sprintf('%d',i)
end
toc
Model_1 = Model_

%% 限制培养条件,测试，仅限制葡萄糖
% 首先关闭大模型中所有大模型中所有反应交换
r_EX_ = find(contains(Model_1.varNames,'EX_'));
Model_1.var_lb(r_EX_,:) = 0;

% 确定碳源
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


% 打开基础物质的交换反应：
wait_open_rxn(4,:) = []
for i = 1:nutrition_conditions
    for ii = 1:length(wait_open_rxn)
        c_ = find(strcmp(Model_1.varNames,strcat(num2str(i),'_',wait_open_rxn{ii,2})));
        Model_1.var_lb(c_,1) = -100;
    end
    clc
    sprintf('已经进行了%d',i)
end

% load need_fill_gap.mat
c_ = find(strcmp(Model_1.varNames,'1_EX_cpd00232_e0'));
Model_1.var_lb(c_,1) = -100;

Model_2 = Model_1
rr = find(contains(Model_2.varNames,'_add_CG_biomass cgl ATCC13032'))
Model_2.f(rr,1) = 1
solution = gurobi_solve_CY(Model_2,'max') 
Model_2.f = zeros(length(Model_2.f),1)

%% 设置目标函数
load recorde_rxn_score.mat
data_rxn.normalization = normalize(data_rxn.bitscore,'range')-1

 
r1 = find(~contains(Universal_Model.rxns,'EX_'));
r2 = find(~contains(Universal_Model.rxns,'biomass'));
r = intersect(r1,r2);
rxn_wait_selected = Universal_Model.rxns(r,:);

% rep和add直接给一个高得分即可 (处理初始模型)

f = @(x)strrep(x,'_c0','')
rxn_wait_selected(:,1) = cellfun(f,rxn_wait_selected(:,1),'UniformOutput',false)
r_rep = find(contains(rxn_wait_selected(:,1),'rep_'));
r_add = find(contains(rxn_wait_selected(:,1),'add_'));
r_ = union(r_rep,r_add);
for i = 1:length(r_)
    rxn_wait_selected{r_(i),2} = 0.1;
end

% 处理自发反应()
r_spontaneous = find(strcmp(data_rxn.status,'spontaneous'))
data_rxn_spontaneous = data_rxn(r_spontaneous,:)

for i = 1:length(rxn_wait_selected)
    if length(rxn_wait_selected{i,2})==0
        if sum(contains(data_rxn_spontaneous.dbhit,rxn_wait_selected{i,1}))>0
            rxn_wait_selected{i,2} = 0;
            sprintf(rxn_wait_selected{i,1});
        end
    end
end     % 发现好像没有自发反应：

% 给反应赋予相应分值：
data_rxn = sortrows(data_rxn,21,'descend')
for i = 1:length(rxn_wait_selected)
    if length(rxn_wait_selected{i,2})==0
        if sum(contains(data_rxn.dbhit,rxn_wait_selected{i,1}))>0
            r_ = find(contains(data_rxn.dbhit,rxn_wait_selected{i,1}));
            rxn_wait_selected{i,2} = data_rxn.normalization(r_(1),1);
            sprintf(rxn_wait_selected{i,1})
        else
            rxn_wait_selected{i,2} = min(data_rxn.normalization);
        end 
    end
end

%% 设置比生长速率下限和约束类型：
for i = 1:nutrition_conditions
    c_ = find(strcmp(strcat(num2str(i),'_add_CG_biomass cgl ATCC13032'),Model_2.varNames));
    Model_2.var_lb(c_,1) = 10;
    clc
    sprintf('已经进行了%d',i)
end

for i = 14369:length(Model_2.varNames)
    str_1 = strrep(Model_2.varNames{i,1},'select_or_not','');
    str_1 = strrep(str_1,'_c0','');
    str_1 = strrep(str_1,'_f','');
    str_1 = strrep(str_1,'_r','');
    r__ = find(strcmp(rxn_wait_selected(:,1),str_1));
    Model_2.f(i,1) = rxn_wait_selected{r__,2};
end


solution = gurobi_solve_CY(Model_2,'max') 

r_ = find(solution.x~=0)
Res = [Model_2.varNames(r_,1),num2cell(solution.x(r_,1))]



%% MILP 矩阵初步构建
% 先尝试1培养基的条件,这次反应得分要用 zscore的方法进行归一化处理，并且拆分方程的可逆性。
load Unversal-model_icw773metadd_rxnadd_rxnrep.mat

nutrition_conditions = 1;
Model.A = sparse(nutrition_conditions.*size(Universal_Model.S,1),nutrition_conditions.*size(Universal_Model.S,2));
Model.constraintType = {};
Model.constraintNames = {};
Model.rhs = [];
b = 1;
% 增加1种培养条件的约束
for i = 1:nutrition_conditions
    for ii=1:size(Universal_Model.mets,1)
        Model.constraintType{b,1} = '=';
        Model.constraintNames{b,1} = strcat(num2str(i),'_',Universal_Model.mets{ii});  % 代谢物 质量平衡约束
        Model.rhs(b,1) = 0;
        Model.A(b,:) = Universal_Model.S(ii,:);
        b = b + 1;
    end
end

% 增加1种培养条件的变量
Model.var_lb = [];
Model.var_ub = [];
Model.varNames = {}
Model.vartypes = []
for i = 1:nutrition_conditions
    for ii=1:size(Universal_Model.rxns,1)
        Model.var_lb = [Model.var_lb;Universal_Model.lb(ii)];
        Model.var_ub = [Model.var_ub;Universal_Model.ub(ii)];
        Model.varNames = [Model.varNames;strcat(num2str(i),'_',Universal_Model.rxns{ii})];
        Model.vartypes = [Model.vartypes;'C'];
    end
    clc
    sprintf('已经进行了%d',i)
end
Model.f = zeros(size(Model.A,2),1);
Model.vartypes = num2cell(Model.vartypes)

% 增加二元变量，用来限制反应数量
r1 = find(~contains(Universal_Model.rxns,'EX_'));
r2 = find(~contains(Universal_Model.rxns,'biomass'));
r = intersect(r1,r2)
rxn_wait_selected = Universal_Model.rxns(r,:);
for i = 1:length(r)
    Model = addNewVariable_CY(Model, strcat('select_or_not',rxn_wait_selected{i},'_f'),'B',[0,1]);
    Model = addNewVariable_CY(Model, strcat('select_or_not',rxn_wait_selected{i},'_r'),'B',[0,1]);
    clc
    sprintf('已经进行了%d',i)
end

tic
for i = 1:nutrition_conditions
    for ii = 1:length(r)
        Model.constraintNames = [Model.constraintNames;strcat('Constraint_lb_',num2str(i),'_',Universal_Model.rxns{r(ii),1})];
        Model.constraintNames = [Model.constraintNames;strcat('Constraint_ub_',num2str(i),'_',Universal_Model.rxns{r(ii),1})];
        Model.constraintNames = [Model.constraintNames;strcat('Constraint_ub_lb',num2str(i),'_',Universal_Model.rxns{r(ii),1})];
        Model.constraintType = [Model.constraintType;{'>'}];
        Model.constraintType = [Model.constraintType;{'<'}];
        Model.constraintType = [Model.constraintType;{'<'}];
        Model.rhs = [Model.rhs;0;0;1];
    end
    clc
    sprintf('%d',i)
end
toc

% 大矩阵中填写相应的系数
Model_ = Model

tic
for i = 1:nutrition_conditions
    for ii = 1:length(r)
        c_ = find(strcmp(Model_.varNames,strcat(num2str(i),'_',Universal_Model.rxns{r(ii),1})));
        pos_f = find(strcmp(Model_.varNames,strcat('select_or_not',rxn_wait_selected{i},'_f')));
        pos_r = find(strcmp(Model_.varNames,strcat('select_or_not',rxn_wait_selected{i},'_r')));
        a = sparse(3,size(Model_.A,2));
        CLHS.pos = [c_,pos_f,pos_r];
        CLHS.coeffs_1 = [1,-1e-6,1e6 ];
        CLHS.coeffs_2 = [1,-1e6,1e-6];
        a(1,CLHS.pos) = CLHS.coeffs_1;
        a(2,CLHS.pos) = CLHS.coeffs_2;
        a(3,[pos_f,pos_r]) = 1;
        Model_.A = [Model_.A;a];
    end
    clc
    sprintf('%d',i)
end
toc
Model_1 = Model_

%% 限制培养条件,测试，仅限制葡萄糖
% 首先关闭大模型中所有大模型中所有反应交换
r_EX_ = find(contains(Model_1.varNames,'EX_'));
Model_1.var_lb(r_EX_,:) = 0;

% 确定碳源
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


% 打开基础物质的交换反应：
wait_open_rxn(4,:) = []
for i = 1:nutrition_conditions
    for ii = 1:length(wait_open_rxn)
        c_ = find(strcmp(Model_1.varNames,strcat(num2str(i),'_',wait_open_rxn{ii,2})));
        Model_1.var_lb(c_,1) = -100;
    end
    clc
    sprintf('已经进行了%d',i)
end

% load need_fill_gap.mat
c_ = find(strcmp(Model_1.varNames,'1_EX_cpd00232_e0'));
Model_1.var_lb(c_,1) = -100;

Model_2 = Model_1
%% 设置目标函数
load recorde_rxn_score.mat
data_rxn.normalization = normalize(data_rxn.bitscore,'range')-1

 
r1 = find(~contains(Universal_Model.rxns,'EX_'));
r2 = find(~contains(Universal_Model.rxns,'biomass'));
r = intersect(r1,r2);
rxn_wait_selected = Universal_Model.rxns(r,:);

% rep和add直接给一个高得分即可 (处理初始模型)

f = @(x)strrep(x,'_c0','')
rxn_wait_selected(:,1) = cellfun(f,rxn_wait_selected(:,1),'UniformOutput',false)
r_rep = find(contains(rxn_wait_selected(:,1),'rep_'));
r_add = find(contains(rxn_wait_selected(:,1),'add_'));
r_ = union(r_rep,r_add);
for i = 1:length(r_)
    rxn_wait_selected{r_(i),2} = 0.1;
end

% 处理自发反应()
r_spontaneous = find(strcmp(data_rxn.status,'spontaneous'))
data_rxn_spontaneous = data_rxn(r_spontaneous,:)

for i = 1:length(rxn_wait_selected)
    if length(rxn_wait_selected{i,2})==0
        if sum(contains(data_rxn_spontaneous.dbhit,rxn_wait_selected{i,1}))>0
            rxn_wait_selected{i,2} = 0;
            sprintf(rxn_wait_selected{i,1});
        end
    end
end     % 发现好像没有自发反应：

% 给反应赋予相应分值：
data_rxn = sortrows(data_rxn,21,'descend')
for i = 1:length(rxn_wait_selected)
    if length(rxn_wait_selected{i,2})==0
        if sum(contains(data_rxn.dbhit,rxn_wait_selected{i,1}))>0
            r_ = find(contains(data_rxn.dbhit,rxn_wait_selected{i,1}));
            rxn_wait_selected{i,2} = data_rxn.normalization(r_(1),1);
            sprintf(rxn_wait_selected{i,1})
        else
            rxn_wait_selected{i,2} = min(data_rxn.normalization);
        end 
    end
end

%% 设置比生长速率下限和约束类型：
for i = 1:nutrition_conditions
    c_ = find(strcmp(strcat(num2str(i),'_add_CG_biomass cgl ATCC13032'),Model_2.varNames));
    Model_2.var_lb(c_,1) = 10;
    clc
    sprintf('已经进行了%d',i)
end

for i = 14369:length(Model_2.varNames)
    str_1 = strrep(Model_2.varNames{i,1},'select_or_not','');
    str_1 = strrep(str_1,'_c0','');
    str_1 = strrep(str_1,'_f','');
    str_1 = strrep(str_1,'_r','');
    r__ = find(strcmp(rxn_wait_selected(:,1),str_1));
    Model_2.f(i,1) = rxn_wait_selected{r__,2};
end


solution = gurobi_solve_CY(Model_2,'max') 

r_ = find(solution.x~=0)
Res = [Model_2.varNames(r_,1),num2cell(solution.x(r_,1))]


%% 限制5种不同的培养条件：
load Unversal-model_icw773metadd_rxnadd_rxnrep.mat

nutrition_conditions = 5;
Model.A = sparse(nutrition_conditions.*size(Universal_Model.S,1),nutrition_conditions.*size(Universal_Model.S,2));
Model.constraintType = {};
Model.constraintNames = {};
Model.rhs = [];
b = 1;
% 增加1种培养条件的约束
for i = 1:nutrition_conditions
    for ii=1:size(Universal_Model.mets,1)
        Model.constraintType{b,1} = '=';
        Model.constraintNames{b,1} = strcat(num2str(i),'_',Universal_Model.mets{ii});  % 代谢物 质量平衡约束
        Model.rhs(b,1) = 0;
        b = b + 1;
    end
end

% 增加1种培养条件的变量
Model.var_lb = [];
Model.var_ub = [];
Model.varNames = {}
Model.vartypes = []
for i = 1:nutrition_conditions
    for ii=1:size(Universal_Model.rxns,1)
        Model.var_lb = [Model.var_lb;Universal_Model.lb(ii)];
        Model.var_ub = [Model.var_ub;Universal_Model.ub(ii)];
        Model.varNames = [Model.varNames;strcat(num2str(i),'_',Universal_Model.rxns{ii})];
        Model.vartypes = [Model.vartypes;'C'];
    end
    clc
    sprintf('已经进行了%d',i)
end
Model.f = zeros(size(Model.A,2),1);
Model.vartypes = num2cell(Model.vartypes)

% 增加二元变量，用来限制反应数量
r1 = find(~contains(Universal_Model.rxns,'EX_'));
r2 = find(~contains(Universal_Model.rxns,'biomass'));
r = intersect(r1,r2)
rxn_wait_selected = Universal_Model.rxns(r,:);
for i = 1:length(r)
    Model = addNewVariable_CY(Model, strcat('select_or_not',rxn_wait_selected{i}),'B',[0,1]);
    clc
    sprintf('已经进行了%d',i)
end

tic
for i = 1:nutrition_conditions
    for ii = 1:length(r)
        Model.constraintNames = [Model.constraintNames;strcat('Constraint_lb_',num2str(i),'_',Universal_Model.rxns{r(ii),1})];
        Model.constraintNames = [Model.constraintNames;strcat('Constraint_ub_',num2str(i),'_',Universal_Model.rxns{r(ii),1})];
        Model.constraintType = [Model.constraintType;{'>'}];
        Model.constraintType = [Model.constraintType;{'<'}];
        Model.rhs = [Model.rhs;0;0];
    end
    clc
    sprintf('%d',i)
end
toc

% 大矩阵中填写相应的系数
Model_ = Model
tic
for i = 1:nutrition_conditions
    for ii = 1:length(r)
        c_ = find(strcmp(Model_.varNames,strcat(num2str(i),'_',Universal_Model.rxns{r(ii),1})));
        pos_b = find(strcmp(Model_.varNames,strcat('select_or_not',rxn_wait_selected{i})));                                                % 修正
        a = sparse(2,size(Model_.A,2));
        CLHS.pos = [c_,pos_b];
        CLHS.coeffs_1 = [1,1e6];
        CLHS.coeffs_2 = [1,-1e6];
        a(1,CLHS.pos) = CLHS.coeffs_1;
        a(2,CLHS.pos) = CLHS.coeffs_2;
        Model_.A = [Model_.A;a];
    end
    clc
    sprintf('%d',i)
end
toc
Model_1 = Model_

%% 限制培养条件,测试，仅限制葡萄糖
% 首先关闭大模型中所有大模型中所有反应交换
r_EX_ = find(contains(Model_1.varNames,'EX_'));
Model_1.var_lb(r_EX_,:) = 0;

% 确定基本物质
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

% 确定相应碳源
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
    r = find(strcmp(carbon_source_substance{i,1},Universal_Model.mets));
    c = find(Universal_Model.S(r,:)==-1);
    for ii = 1:length(c)
        if length(find(Universal_Model.S(:,c(ii))~=0))==1;
            carbon_source_substance{i,4} = Universal_Model.rxns{c(ii),1};
        end
    end
end

r = find(strcmp(carbon_source_substance(:,3),'False')==1)
carbon_source_substance(r,:) = []
carbon_source_substance = carbon_source_substance(1:5,:)

% 打开基础物质的交换反应：
for i = 1:nutrition_conditions
    c_0 = find(strcmp(Model_1.varNames,strcat(num2str(i),'_',carbon_source_substance{i,4})));
    Model_1.var_lb(c_0,1) = -100;
    for ii = 1:length(wait_open_rxn)
        c_ = find(strcmp(Model_1.varNames,strcat(num2str(i),'_',wait_open_rxn{ii,2})));
        Model_1.var_lb(c_,1) = -100;
    end
    clc
    sprintf('已经进行了%d',i)
end
Model_2 = Model_1

%% 设置目标函数
load recorde_rxn_score.mat
data_rxn.normalization = normalize(data_rxn.bitscore,'range')-1


r1 = find(~contains(Universal_Model.rxns,'EX_'));
r2 = find(~contains(Universal_Model.rxns,'biomass'));
r = intersect(r1,r2);
rxn_wait_selected = Universal_Model.rxns(r,:);

% rep和add直接给一个高得分即可 (处理初始模型)

f = @(x)strrep(x,'_c0','')
rxn_wait_selected(:,1) = cellfun(f,rxn_wait_selected(:,1),'UniformOutput',false)
r_rep = find(contains(rxn_wait_selected(:,1),'rep_'));
r_add = find(contains(rxn_wait_selected(:,1),'add_'));
r_ = union(r_rep,r_add);
for i = 1:length(r_)
    rxn_wait_selected{r_(i),2} = 0.1;
end

% 处理自发反应()
r_spontaneous = find(strcmp(data_rxn.status,'spontaneous'))
data_rxn_spontaneous = data_rxn(r_spontaneous,:)

for i = 1:length(rxn_wait_selected)
    if length(rxn_wait_selected{i,2})==0
        if sum(contains(data_rxn_spontaneous.dbhit,rxn_wait_selected{i,1}))>0
            rxn_wait_selected{i,2} = 0;
            sprintf(rxn_wait_selected{i,1});
        end
    end
end     % 发现好像没有自发反应：

% 给反应赋予相应分值：
data_rxn = sortrows(data_rxn,21,'descend')
for i = 1:length(rxn_wait_selected)
    if length(rxn_wait_selected{i,2})==0
        if sum(contains(data_rxn.dbhit,rxn_wait_selected{i,1}))>0
            r_ = find(contains(data_rxn.dbhit,rxn_wait_selected{i,1}));
            rxn_wait_selected{i,2} = data_rxn.normalization(r_(1),1);
            sprintf(rxn_wait_selected{i,1})
        else
            rxn_wait_selected{i,2} = min(data_rxn.normalization);
        end 
    end
end

%% 设置比生长速率下限和约束类型：
for i = 1:nutrition_conditions
    c_ = find(strcmp(strcat(num2str(i),'_add_CG_biomass cgl ATCC13032'),Model_2.varNames));
    Model_2.var_lb(c_,1) = 10;
    clc
    sprintf('已经进行了%d',i)
end

for i = 71841:length(Model_2.varNames)
    str_1 = strrep(Model_2.varNames{i,1},'select_or_not','');
    str_1 = strrep(str_1,'_c0','');
    r__ = find(strcmp(rxn_wait_selected(:,1),str_1));
    Model_2.f(i,1) = rxn_wait_selected{r__,2};
end

solution = gurobi_solve_CY(Model_2,'max');







%% nchoosek
c = nchoosek([1:20],2)
for i = 1:length(c)
    if sum(solution.pool(c(i,1)).xn == solution.pool(c(i,2)).xn) < length(solution.x)
        c(i,3) = 1;
    else
        c(i,3) = 0;
    end
end                   % 验证完毕，解确实都是完全不相等的。













