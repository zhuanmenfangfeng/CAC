% 经典代谢网络模型集FVA求解，葡萄糖和氧气的吸收速率都设置为10
clear;
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\5.经典代谢网络模型集和多尺度模型集模拟的结果的区别
load oxygen_glu.mat

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

load 500_ensemble.mat

f1 = @(x)x(1:5)
f2 = @(x)strrep(x,'select_or_not_','')

% 设置基础物质的吸收能力
r_EX_ = find(contains(Model.rxns,'EX_'));
Model.lb(r_EX_,:) = 0;

cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\1.代谢网络构建\代谢网络draft
[~,~,basic_substance] = xlsread('biolog_base_composition.csv')
basic_substance = basic_substance(:,2)
basic_substance(1,:) = []
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建

% 由于icw773生物量反应的特殊性，需要额外添加这些物质
waitaddmets_new = {'cpd00244[e0]';'cpd00104[e0]';'cpd00099[e0]';'cpd00027[e0]'}
wait_open_rxn = [waitaddmets_new;basic_substance]
f = @(x)strrep(x,'_e','[e0]');
wait_open_rxn = cellfun(f,wait_open_rxn,'UniformOutput',false);

% 确定对应的交换反应
Universal_Model = Model
for i = 1:length(wait_open_rxn)
    r = find(strcmp(wait_open_rxn{i,1},Universal_Model.mets));
    c = find(Universal_Model.S(r,:)==-1);
    for ii = 1:length(c)
        if length(find(Universal_Model.S(:,c(ii))~=0))==1;
            wait_open_rxn{i,2} = Universal_Model.rxns{c(ii),1};
        end
    end
end

wait_open_rxn(4,:) = [];

% 设置相应的反应速度
for i = 1:length(wait_open_rxn)
    c_ = find(strcmp(Model.rxns,wait_open_rxn{i,2}));
    Model.lb(c_,1) = -100;
end

Model = changeObjective(Model,'add_CG_biomass cgl ATCC13032');
% optimizeCbModel(Model)       小验证一下，果然没有解(因为没有对应的碳源)

load big_matrix_example.mat
Res = Model_.varNames;
Res(:,2) = cellfun(f1,Res(:,1),'UniformOutput',false);
% 待选取的变量
r1 = find(strcmp(Res(:,2),'selec'));

% 登录碳源的培养条件
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\1.代谢网络构建\代谢网络draft;
load finall_gap_fill.mat;
carbon_source_substance = carbon_source_substance(2:6,:);

classical_fva = {};
for i = 1:50
    clearvars -except classical_fva Res rxn_wait_selected x y solution i r1 f2 Model
    Res_1 = Res;
    Res_1(:,3) = num2cell(solution.pool(i).xn);
    r2 = find(cell2mat(Res_1(:,3))==1);
    r = intersect(r1,r2);

    must_need_rxns = Res_1(r,1);
    must_need_rxns = cellfun(f2,must_need_rxns,'UniformOutput',false);
    wait_remove = setdiff(rxn_wait_selected,must_need_rxns);
    Model_1 = removeRxns(Model,wait_remove);
    rxns_list = Model_1.rxns;

    r_oxygen = find(strcmp(Model_1.rxns,'rep_EX_o2(e)'));
    r_glucose = find(strcmp(Model_1.rxns,'rep_EX_glc(e)'));

    Model_1 = changeRxnBounds(Model_1, 'rep_EX_o2(e)', -x(i), 'l');
    Model_1 = changeRxnBounds(Model_1, 'rep_EX_glc(e)', -y(i), 'l');

    for ii = 1:size(rxns_list,1)
        Model_1 = changeObjective(Model_1,rxns_list{ii,1});
        solution_max = optimizeCbModel(Model_1,'max');
        solution_min = optimizeCbModel(Model_1,'min');
        rxns_list{ii,2} = solution_max.f;
        rxns_list{ii,3} = solution_min.f;
    end

    classical_fva{i,1} = rxns_list;
    clc
    sprintf('%d',i)
end

% 最后可以统一插值
fva_mat_ = []
xq = 0:0.0002:1;
for i = 1:20
    data = classical_fva{i,1};
    data_flux = cell2mat(data(:,2:3));
    fva_mat = data_flux(:,1) - data_flux(:,2);
    fva_mat = sort(abs(fva_mat));
    for ii = 1:size(fva_mat,1)
    % 计算下频率
        fva_mat(ii,2) = ii./size(fva_mat,1);
    end
    % 对fva值进行插值
    vq = interp1(fva_mat(:,2),fva_mat(:,1),xq);
    fva_mat_ = [fva_mat_,vq'];
end

% 再整合一个x列
fva_mat_ = [xq',fva_mat_];
fva_mat_(1:4,:) = [];

% 求平均值：
mean_ = mean(fva_mat_(:,2:21),2)
std_ = std(fva_mat_(:,2:21),0,2)

cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\6.模型集FVA求解


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 多尺度模型集模拟 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\5.经典代谢网络模型集和多尺度模型集模拟的结果的区别
load oxygen_glu.mat

clear;
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

cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\5.经典代谢网络模型集和多尺度模型集模拟的结果的区别
load 500_ensemble.mat

f1 = @(x)x(1:5)
f2 = @(x)strrep(x,'select_or_not_','')

% 设置基础物质的吸收能力
r_EX_ = find(contains(Model.rxns,'EX_'));
Model.lb(r_EX_,:) = 0;

cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\1.代谢网络构建\代谢网络draft
[~,~,basic_substance] = xlsread('biolog_base_composition.csv')
basic_substance = basic_substance(:,2)
basic_substance(1,:) = []
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建

% 由于icw773生物量反应的特殊性，需要额外添加这些物质
waitaddmets_new = {'cpd00244[e0]';'cpd00104[e0]';'cpd00099[e0]';'cpd00027[e0]'}
wait_open_rxn = [waitaddmets_new;basic_substance]
f = @(x)strrep(x,'_e','[e0]');
wait_open_rxn = cellfun(f,wait_open_rxn,'UniformOutput',false);

% 确定对应的交换反应
Universal_Model = Model
for i = 1:length(wait_open_rxn)
    r = find(strcmp(wait_open_rxn{i,1},Universal_Model.mets));
    c = find(Universal_Model.S(r,:)==-1);
    for ii = 1:length(c)
        if length(find(Universal_Model.S(:,c(ii))~=0))==1;
            wait_open_rxn{i,2} = Universal_Model.rxns{c(ii),1};
        end
    end
end

wait_open_rxn(4,:) = [];

% 修改相应反应的格式

for i = 1:length(wait_open_rxn)
    wait_open_rxn{i,3} = strrep(wait_open_rxn{i,2},'(e)','_e');
end

cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\5.经典代谢网络模型集和多尺度模型集模拟的结果的区别
load oxygen_glu.mat

Res_f = {}

for i = 1:20
    cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建\Multiscale_Model_ensemble
    load(['model-',num2str(i)]);
    model = Model_thermo_;

    % 首先确定要引入多少变量：
    r1 = find(contains(model.varNames,'F_'));
    r2 = find(~contains(model.varNames,'No'));
    r3 = find(~contains(model.varNames,'FU'));
    r4 = find(~contains(model.varNames,'BU'));
    r5 = find(~contains(model.varNames,'DG'));
    r6 = find(~contains(model.varNames,'DGo'));

    r = intersect(r1,r2);
    r = intersect(r,r3);
    r = intersect(r,r4);
    r = intersect(r,r5);
    r = intersect(r,r6);
    rxns_variable = model.varNames(r,:);

    % 挨个为其创建约束条件
    for ii = 1:length(rxns_variable)
        if sum(strcmp(model.varNames,rxns_variable{ii,1}))>0 & sum(strcmp(model.varNames,['R_',rxns_variable{ii,1}(3:end)]))>0

            c1 = find(strcmp(model.varNames,rxns_variable{ii,1}));
            c2 = find(strcmp(model.varNames,['R_',rxns_variable{ii,1}(3:end)]));
            c = [c1,c2];

            num_vars = length(model.varNames);
            num_constr = length(model.constraintNames);

            model.varNames{num_vars+1,1} = ['Sum_',rxns_variable{ii,1}(3:end)];
            model.f   = [model.f;0];
            model.var_lb(num_vars+1,1)   = -1000;                          % lower bound
            model.var_ub(num_vars+1,1)   = 1000;                       % upper bound10
            model.vartypes = [model.vartypes;'C'];
            model.A(:,num_vars+1) = zeros(num_constr,1);                    % 增加一列变量以表示新合并的变量。
            model.constraintNames{num_constr+1,1} = ['flux_sum',rxns_variable{ii,1}(3:end)];                 % strcat('DFSEU_',model.mets{i});
            model.constraintType = [model.constraintType;'='];
            model.rhs = [model.rhs;0];
            c = [c,num_vars+1];
            model.A(num_constr+1,c) = [-ones(1,length(c)-1),1];   % 分别在约束中填补对应的系数 这个修改好了可以回头设置成函数。
        end
    end

    cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\6.模型集FVA求解
    name = ['FVAmodel-',num2str(i),'.mat'];
    save(name,'model');
end

clear

%% 计算一下FVA:

load FVAmodel-1.mat
gurobi_solve_CY(model,'max')                % 还能求解呢，说明增加的过程对模型没有什么影响。


% 统计一下有多少反应：
r = find(contains(model.varNames,'Sum_'))
rxns_name = model.varNames(r,:)
for i = 1:size(rxns_name,1)
    c_ = find(strcmp(model.varNames,rxns_name{i,1}));
    model.f = zeros(length(model.f),1);
    model.f(c_,1) = 1;
    solution = gurobi_solve_CY(model,'max');
    rxns_name{i,2} = solution.objval;
    solution = gurobi_solve_CY(model,'min');
    rxns_name{i,3} = solution.objval;
    clc
    sprintf('已经进行了%d',i)
end                                         %% 完美运行

%% 设定一下碳源吧：

clear;
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\5.经典代谢网络模型集和多尺度模型集模拟的结果的区别
load oxygen_glu.mat

cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\6.模型集FVA求解

Res_ = {}
for i = 1:20
    load(['FVAmodel-',num2str(i),'.mat']);
    r = find(contains(model.varNames,'Sum_'));
    rxns_name = model.varNames(r,:);
    r_oxygen = find(strcmp(model.varNames,'R_rep_EX_o2_e'));
    r_glucose = find(strcmp(model.varNames,'R_rep_EX_glc_e'));
    model.var_ub(r_oxygen,1) = x(i)
    model.var_ub(r_glucose,1) = y(i)

    for ii = 1:size(rxns_name,1)
        c_ = find(strcmp(model.varNames,rxns_name{ii,1}));
        model.f = zeros(length(model.f),1);
        model.f(c_,1) = 1;
        solution = gurobi_solve_CY(model,'max');
        rxns_name{ii,2} = solution.objval;
        solution = gurobi_solve_CY(model,'min');
        rxns_name{ii,3} = solution.objval;
        clc
        sprintf('已经进行了%d',ii)
    end
    Res_{i,1} = rxns_name;
    clc
    sprintf('已经进行了%d',i)
end

% 先用6个做一下吧：

fva_mat_ = [];
xq = 0:0.0002:1;
for i = 1:6
    data = Res_{i,1};
    data_flux = cell2mat(data(:,2:3));
    fva_mat = data_flux(:,1) - data_flux(:,2);
    fva_mat = sort(abs(fva_mat));
    for ii = 1:size(fva_mat,1)
    % 计算下频率
        fva_mat(ii,2) = ii./size(fva_mat,1);
    end
    % 对fva值进行插值
    vq = interp1(fva_mat(:,2),fva_mat(:,1),xq);
    fva_mat_ = [fva_mat_,vq'];
end

fva_mat_ = [xq',fva_mat_];
fva_mat_(1:4,:) = [];

MEAN_ = mean(fva_mat_(:,2:7),2)
STD_ = std(fva_mat_(:,2:7),0,2)
