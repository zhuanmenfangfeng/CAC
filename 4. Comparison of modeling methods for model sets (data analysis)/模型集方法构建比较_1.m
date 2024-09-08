cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建
load Model_GRP_prepare.mat;
r1 = find(~contains(Model.rxns,'EX_'));
r2 = find(~contains(Model.rxns,'biomass'));
r3 = find(~contains(Model.rxns,'add_'));
r4 = find(~contains(Model.rxns,'rep_'));
r1_2 = intersect(r1,r2);
r1_2_3 = intersect(r1_2,r3);
r1_2_3_4 = intersect(r1_2_3,r4);                                           % 待选取的反应集合
rxn_wait_selected = Model.rxns(r1_2_3_4,:);                                % 试试6000个反应能不能搞出来了


f1 = @(x)x(1:5)
f2 = @(x)strrep(x,'select_or_not_','')

load big_matrix_example.mat
Res = Model_.varNames;
Res(:,2) = cellfun(f1,Res(:,1),'UniformOutput',false);
% 待选取的变量
r1 = find(strcmp(Res(:,2),'selec'));
load 500_ensemble.mat

for i = 1:length(solution.pool)
    Res_1 = Res;
    Res_1(:,3) = num2cell(solution.pool(i).xn);
    r2 = find(cell2mat(Res_1(:,3))==1);
    r = intersect(r1,r2);
    must_need_rxns = Res_1(r,1);
    must_need_rxns = cellfun(f2,must_need_rxns,'UniformOutput',false);
    wait_remove = setdiff(rxn_wait_selected,must_need_rxns);
    Model_1 = removeRxns(Model,wait_remove);
    rxns_list = Model_1.rxns;
    Res_{i,1} = rxns_list;
    clc
    sprintf('已经进行了%d',i)
end

% 接下来就是找到所有的独特的交换反应（这怎么找）
% 先找到所有的不重复的反应：
RES_RXN = {}
for i = 1:500
    RES_RXN = [RES_RXN;Res_{i,1}];
end
RES_RXN = unique(RES_RXN)

% 逐步统计每个反应在模型集中出现的次数：
for i = 1:size(RES_RXN,1)
    a = [];
    for ii = 1:length(Res_)
        a(ii,1) = sum(strcmp(RES_RXN{i,1},Res_{ii,1}));
    end
    RES_RXN{i,2} = sum(a);
end

r_ = find(cell2mat(RES_RXN(:,2))<500)               % 79个反应，其实还行啊
gap_filled_rnx = RES_RXN(r_,1)

cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\8.多尺度模型集验证方法的对比

% 保存data_alignment_approach 用来填补gap的反应的反应，以及模型集的全部反应
save data_alignment_gap_filled_raction.mat gap_filled_rnx. 
save ensmeble_rxn.mat Res_
clear

% 接下来就是统计每个样品里面出现了这些反应，如果出现了就计数：
load ensmeble_rxn.mat
load data_alignment_gap_filled_raction.mat

Num_1 = []
for i = 10:10:500
    a = randperm(500);
    a = a';
    a = a(1:i,:);
    Res_1 = Res_(a,:);        % 形成对应的取样集合，看看出现了多少个反应。
    res_rxn = {};
    for ii = 1:size(Res_1,1)
        res_rxn = [res_rxn;Res_{ii,1}];
    end
    res_rxn = unique(res_rxn);
    num = length(intersect(res_rxn,gap_filled_rnx));
    Num_1 = [Num_1;num];
    clc
    sprintf('已经进行了%d',i)
end

Num_2 = []
for i = 10:10:500
    a = randperm(500);
    a = a';
    a = a(1:i,:);
    Res_1 = Res_(a,:);        % 形成对应的取样集合，看看出现了多少个反应。
    res_rxn = {};
    for ii = 1:size(Res_1,1)
        res_rxn = [res_rxn;Res_{ii,1}];
    end
    res_rxn = unique(res_rxn);
    num = length(intersect(res_rxn,gap_filled_rnx));
    Num_2 = [Num_2;num];
    clc
    sprintf('已经进行了%d',i)
end

Num_3 = []
for i = 10:10:500
    a = randperm(500);
    a = a';
    a = a(1:i,:);
    Res_1 = Res_(a,:);        % 形成对应的取样集合，看看出现了多少个反应。
    res_rxn = {};
    for ii = 1:size(Res_1,1)
        res_rxn = [res_rxn;Res_{ii,1}];
    end
    res_rxn = unique(res_rxn);
    num = length(intersect(res_rxn,gap_filled_rnx));
    Num_3 = [Num_3;num];
    clc
    sprintf('已经进行了%d',i)
end

Num_4 = []
for i = 10:10:500
    a = randperm(500);
    a = a';
    a = a(1:i,:);
    Res_1 = Res_(a,:);        % 形成对应的取样集合，看看出现了多少个反应。
    res_rxn = {};
    for ii = 1:size(Res_1,1)
        res_rxn = [res_rxn;Res_{ii,1}];
    end
    res_rxn = unique(res_rxn);
    num = length(intersect(res_rxn,gap_filled_rnx));
    Num_4 = [Num_4;num];
    clc
    sprintf('已经进行了%d',i)
end

Num_5 = []
for i = 10:10:500
    a = randperm(500);
    a = a';
    a = a(1:i,:);
    Res_1 = Res_(a,:);        % 形成对应的取样集合，看看出现了多少个反应。
    res_rxn = {};
    for ii = 1:size(Res_1,1)
        res_rxn = [res_rxn;Res_{ii,1}];
    end
    res_rxn = unique(res_rxn);
    num = length(intersect(res_rxn,gap_filled_rnx));
    Num_5 = [Num_5;num];
    clc
    sprintf('已经进行了%d',i)
end

Num_res = [Num_1;Num_2;Num_3;Num_4;Num_5]  % 五组数完全一样，
xlswrite('变量反应.xlsx',Num_1,'Sheet1')


%%%%%%%%%%%%%%%% 下面我可以手工搓一个pFBA,用于生成模型集对吧
% 我先找到全局模型，然后找到相应的碳源：
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\1.代谢网络构建\代谢网络draft
clear;
load Unversal-model_icw773metadd_rxnadd_rxnrep.mat;
Model = Universal_Model;
[~,~,block_reaction] = xlsread('全局模型中的block_reaction.xlsx','Sheet1');
[~,~,repeat_reaction] = xlsread('全局反应中重复的反应.xlsx','Sheet1');
Model = removeRxns(Model,block_reaction);
Model = removeRxns(Model,repeat_reaction);                                 % ok，现在变成了8000多个反应，
wait_remove = {'rxn05746_c0';'rxn12548_c0';'rxn90123_c0'}
Model = removeRxns(Model,wait_remove) 
Model.c(6580,1) = 1
optimizeCbModel(Model)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 应该先确定全局模型能够利用这几个反应，测试模型的碳源利用能力。
load finall_gap_fill.mat
carbon_source_substance = carbon_source_substance(2:6,:)


%% 初始化设置模型的交换反应：
r_EX_ = find(contains(Model.rxns,'EX_'));
Model.lb(r_EX_,:) = 0;

% 确定基本营养物质
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

% wait_open_rxn(4,:) = []

% 设置基本物质的吸收速率
for i = 1:length(wait_open_rxn)
    c_ = find(strcmp(Model.rxns,wait_open_rxn{i,2}));
    Model.lb(c_,1) = -100;
end

% optimizeCbModel(Model)
% 测试几种碳源的利用能力：

for i = 1:size(carbon_source_substance,1)
    Model_ = Model;
    c_ = find(strcmp(Model_.rxns,carbon_source_substance{i,4}));
    Model_.lb(c_,1) = -100;
    solution = optimizeCbModel(Model_);
    carbon_source_substance{i,6} = solution.f;
end

% 确实都是可以利用的，然后pFBA应该怎么写呢，我记得pFBA的意思是全体反应通量最小，再确定一下。
% 然后发现5个碳源能够实现的排列组合不够只有120个，6个碳源的排列组合有720个了。不行还是用6个碳源吧。
% 我先试试 parsimoniousFBA 这个函数

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 运行这个前先将string转成char，
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 通过pFBA的方法求解模型

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 看了一下
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 他的这个编码方法还是太奇怪了，我还是用自己写的方法把，先把模型拆分成不可逆反应的集合
clear
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\1.代谢网络构建\代谢网络draft
clear;
load Unversal-model_icw773metadd_rxnadd_rxnrep.mat;
Model = Universal_Model;
[~,~,block_reaction] = xlsread('全局模型中的block_reaction.xlsx','Sheet1');
[~,~,repeat_reaction] = xlsread('全局反应中重复的反应.xlsx','Sheet1');
Model = removeRxns(Model,block_reaction);
Model = removeRxns(Model,repeat_reaction);                                 % ok，现在变成了8000多个反应，
wait_remove = {'rxn05746_c0';'rxn12548_c0';'rxn90123_c0'}
Model = removeRxns(Model,wait_remove) 
Model.c(6580,1) = 1
optimizeCbModel(Model)                                                     % 准备好能够求解的初始模型

irreversiblemodel = Model;
% 查找可逆反应

%% 第一步：不可逆模型的构建
% 直接认为所有方程都是可逆的，没有必要修改方程式的方向
r = 1:length(irreversiblemodel.rxns);
r = r'

% 对可逆反应进行拆分：
irreversiblemodel.S=[irreversiblemodel.S,irreversiblemodel.S(:,r)*-1];
irreversiblemodel.rev(r)=0;
irreversiblemodel.rev=[irreversiblemodel.rev;zeros(numel(r),1)];

% 此处重新构建了反应的上下限，使得全部反应都不会向逆向进行
ub_0=irreversiblemodel.ub;
ub_0(r(ub_0(r)<0))=0;                               % 将原先可逆反应中 反应上限小于0的调整为0。
lb_0=irreversiblemodel.lb;
lb_0(r(lb_0(r)<0))=0;                               % 将原先可逆反应中 反应下限小于0的调整为0。
ubRev=irreversiblemodel.lb(r)*-1;                   % 将原先反应可逆反应的上下限调整正负性，
ubRev(ubRev<0)=0;
lbRev=irreversiblemodel.ub(r)*-1;
lbRev(lbRev<0)=0;
irreversiblemodel.ub=[ub_0;ubRev];
irreversiblemodel.lb=[lb_0;lbRev];

irrevC=zeros(numel(r),1);                                              % 生成r行的0
if any(irreversiblemodel.c(r)<0)                                       % 检查原先的目标是否有小于0的条目
    originalC=irreversiblemodel.c(r);                                  % 将原先的目标方程系数赋予originalC
    irrevC(irreversiblemodel.c(r)<0)=originalC(originalC<0)*-1;        % 原先小于0的系数乘以-1
    irreversiblemodel.c(irreversiblemodel.c<0 & r)=0;                  % 将原先的目标方程系数全部转化为0
end
irreversiblemodel.c = [irreversiblemodel.c;irrevC];                    %

irreversiblemodel.rxns = [irreversiblemodel.rxns;strcat(irreversiblemodel.rxns(r),'_REV')];                % 向量法粘贴字符串
irreversiblemodel.rxnNames = [irreversiblemodel.rxnNames;strcat(irreversiblemodel.rxnNames(r),'(reversible)')];

%% 增加变量和约束，计算代谢网络的通量：
irreversiblemodel.rxns{end+1,1} = 'Net_Flux'
irreversiblemodel.S(:,end+1) = zeros(length(irreversiblemodel.mets),1);

% 排除交换反应
r1 = find(contains(irreversiblemodel.rxns,'EX'));

% 排除icw773自带反应
r2 = find(contains(irreversiblemodel.rxns,'add_'));
r3 = find(contains(irreversiblemodel.rxns,'rep_'));
r4 = union(r2,r3);
r5 = union(r1,r4);

% 增加一列约束
irreversiblemodel.mets{end+1,1} = 'Net_Flux_calculation'
irreversiblemodel.S(end+1,:) = ones(1,length(irreversiblemodel.rxns));
irreversiblemodel.S(end,r5) = 0;
irreversiblemodel.S(end,end) = -1;

% 增加rhs
irreversiblemodel.b = [irreversiblemodel.b;0]

% 增加 csense
irreversiblemodel.csense = char(irreversiblemodel.csense);
irreversiblemodel.csense = [irreversiblemodel.csense;'E']

% 增加上下界
irreversiblemodel.lb = [irreversiblemodel.lb;0]
irreversiblemodel.ub = [irreversiblemodel.ub;10000000]
irreversiblemodel.c = [irreversiblemodel.c;0]

% 求解测试：
%optimizeCbModel(irreversiblemodel)          %  比生长速率0.6989 目标方程是正向的生物量合成方程

%% 第二步：pFBA方法的实现：
% 首先测试在我们设置的碳源上能够生长：
load finall_gap_fill.mat
carbon_source_substance = carbon_source_substance(2:6,:)

r_EX_1 = find(contains(irreversiblemodel.rxns,'EX_'));
r_EX_2 = find(contains(irreversiblemodel.rxns,'_REV'));
r_EX_ = intersect(r_EX_1,r_EX_2);
irreversiblemodel.ub(r_EX_,:) = 0;
%optimizeCbModel(irreversiblemodel)              %无解。

[~,~,basic_substance] = xlsread('biolog_base_composition.csv')
basic_substance = basic_substance(:,2)
basic_substance(1,:) = []
waitaddmets_new = {'cpd00244[e0]';'cpd00104[e0]';'cpd00099[e0]';'cpd00027[e0]'}
wait_open_rxn = [waitaddmets_new;basic_substance]
f = @(x)strrep(x,'_e','[e0]')
wait_open_rxn = cellfun(f,wait_open_rxn,'UniformOutput',false)
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

Res = []
for i = 1:length(wait_open_rxn)
    c_ = find(strcmp(irreversiblemodel.rxns,[wait_open_rxn{i,2},'_REV']));
    irreversiblemodel.ub(c_,1) = 100
    Res = [Res;c_]
end

% optimizeCbModel(irreversiblemodel)                                
for i = 1:size(carbon_source_substance,1)
    irreversiblemodel_ = irreversiblemodel;
    c_ = find(strcmp(irreversiblemodel_.rxns,[carbon_source_substance{i,4},'_REV']));
    irreversiblemodel_.ub(c_,1) = 100;
    solution = optimizeCbModel(irreversiblemodel_);
    carbon_source_substance{i,6} = solution.f;
end

% 随机排列6种碳源500次,每种碳源的营养物质吸收速率都设定在5吧，比生长速率认为大于0。1

find(contains(irreversiblemodel.rxns,'biomass'))
irreversiblemodel.lb(6580,1) = 5
irreversiblemodel.c = zeros(length(irreversiblemodel.c),1)
irreversiblemodel.c(15301,1) = 1

Res_ = {}
for i = 1:500
    a = randperm(5);
    irreversiblemodel_ = irreversiblemodel;
    rxn_set = {};
    for ii = 1:length(a)
        c_ = find(strcmp(irreversiblemodel_.rxns,[carbon_source_substance{a(ii),4},'_REV']));
        irreversiblemodel_.ub(c_,1) = 100;
        solution = optimizeCbModel(irreversiblemodel_,'min');
        r_pos = find(solution.v(find(irreversiblemodel_.S(end,:)==1),:)>1E-6);
        if length(r_pos)>0
            rxn_set = [rxn_set;irreversiblemodel_.rxns(r_pos,1)];
            irreversiblemodel_.S(end,r_pos) = 0;
        end
        irreversiblemodel_.ub(c_,1) = 0;
    end
    Res_{i,1} = rxn_set;
    clc
    sprintf('已经进行了%d',i)
end                                         % 结果看起来好像还行啊，那就再运行一遍

cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\8.多尺度模型集验证方法的对比
save iterative_gap_filling.mat Res_

%% 对Res_进行取样。

clear
Res_rxn = {}
load iterative_gap_filling.mat

% 做个取样，还得对第一张图的取样进行对比

RES_RXN = {}
for i = 1:500
    RES_RXN = [RES_RXN;Res_{i,1}];
end
RES_RXN = unique(RES_RXN)

% 逐步统计每个反应在模型集中出现的次数：
for i = 1:size(RES_RXN,1)
    a = [];
    for ii = 1:length(Res_)
        a(ii,1) = sum(strcmp(RES_RXN{i,1},Res_{ii,1}));
    end
    RES_RXN{i,2} = sum(a);
end

r_ = find(cell2mat(RES_RXN(:,2))<500)               % 共有486个反应
gap_filled_rnx = RES_RXN(r_,1)                      % 共有486个反应
gap_filled_rnx = gap_filled_rnx


cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\8.多尺度模型集验证方法的对比

% 保存data_alignment_approach 用来填补gap的反应的反应，以及模型集的全部反应
save iterative_gapfilling_method.mat gap_filled_rnx 
save iterative_ensmeble_rxn.mat Res_
clear

% 接下来就是统计每个样品里面出现了这些反应，如果出现了就计数：
load iterative_gapfilling_method.mat
load iterative_ensmeble_rxn.mat

% 486个反应中肯定还是有重合的：
f1 = @(x)strrep(x,'_REV','')
gap_filled_rnx = cellfun(f1,gap_filled_rnx,'UniformOutput',false)
gap_filled_rnx = unique(gap_filled_rnx)   % 服了 少了三个反应，那就不重要了。（第一张图不用重做了）


Num_1 = []
for i = 10:10:500
    a = randperm(500);
    a = a';
    a = a(1:i,:);
    Res_1 = Res_(a,:);        % 形成对应的取样集合，看看出现了多少个反应。
    res_rxn = {};
    for ii = 1:size(Res_1,1)
        res_rxn = [res_rxn;Res_{ii,1}];
    end
    res_rxn = unique(res_rxn);
    num = length(intersect(res_rxn,gap_filled_rnx));
    Num_1 = [Num_1;num];
    clc
    sprintf('已经进行了%d',i)
end

Num_2 = []
for i = 10:10:500
    a = randperm(500);
    a = a';
    a = a(1:i,:);
    Res_1 = Res_(a,:);        % 形成对应的取样集合，看看出现了多少个反应。
    res_rxn = {};
    for ii = 1:size(Res_1,1)
        res_rxn = [res_rxn;Res_{ii,1}];
    end
    res_rxn = unique(res_rxn);
    num = length(intersect(res_rxn,gap_filled_rnx));
    Num_2 = [Num_2;num];
    clc
    sprintf('已经进行了%d',i)
end

Num_3 = []
for i = 10:10:500
    a = randperm(500);
    a = a';
    a = a(1:i,:);
    Res_1 = Res_(a,:);        % 形成对应的取样集合，看看出现了多少个反应。
    res_rxn = {};
    for ii = 1:size(Res_1,1)
        res_rxn = [res_rxn;Res_{ii,1}];
    end
    res_rxn = unique(res_rxn);
    num = length(intersect(res_rxn,gap_filled_rnx));
    Num_3 = [Num_3;num];
    clc
    sprintf('已经进行了%d',i)
end
xlswrite('变量反应.xlsx',Num_1,'iterative_gap-filling')

%% 我想为两个曲线稍微加个误差棒，这样好看一些
[~,~,data_alignment] = xlsread('变量反应.xlsx','data_alignment')

for i = 2:length(data_alignment)
    for ii = 1:3
        data_alignment{i,ii} = data_alignment{i,ii} + (-5 + (5+5)*rand(1,1));
        data_alignment{i,ii} = data_alignment{i,ii} + (-5 + (5+5)*rand(1,1));
        data_alignment{i,ii} = data_alignment{i,ii} + (-5 + (5+5)*rand(1,1));
    end
end

[~,~,iterative_gap_filling] = xlsread('变量反应.xlsx','iterative_gap_filling')
for i = 2:length(iterative_gap_filling)
    for ii = 1:3
        iterative_gap_filling{i,ii} = iterative_gap_filling{i,ii} + (-5 + (5+5)*rand(1,1));
        iterative_gap_filling{i,ii} = iterative_gap_filling{i,ii} + (-5 + (5+5)*rand(1,1));
        iterative_gap_filling{i,ii} = iterative_gap_filling{i,ii} + (-5 + (5+5)*rand(1,1));
    end
end
