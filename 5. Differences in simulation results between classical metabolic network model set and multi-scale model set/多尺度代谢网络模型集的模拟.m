% 验证一下对不同培养条件的利用能力：

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

Res_f = {}
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建\Multiscale_Model_ensemble
for i = 1:500
    load(['model-',num2str(i)]);
    model = Model_thermo_;
    r_EX_ = find(contains(model.varNames,'EX_'));
    model.var_lb(r_EX_,:) = 0;
    for ii = 1:length(wait_open_rxn)
        c_ = find(strcmp(model.varNames,['R_',wait_open_rxn{ii,3}]));
        model.var_lb(c_,1) = -100;
    end
    r_oxygen = find(strcmp(model.varNames,'R_rep_EX_o2_e'));
    r_glucose = find(strcmp(model.varNames,'R_rep_EX_glc_e'));
    b = 1;
    for ii = 0:3:12
        for iii = 0:3:12
            model.var_ub(r_oxygen,1) = ii;
            model.var_ub(r_glucose,1) = iii;
            solution_1 = gurobi_solve_CY(model,'max')
            Res_f{i,b} = solution_1.objval;
            b = b + 1;
        end
    end
    clc
    sprintf('正在处理模型%d',i)
end
Res_f_m = Res_f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 聚类结果分析 %%%%%%%%%%%%%%%%%%%%%%%%%%

data = [Res_f;Res_f_m]
data = cell2mat(data)
[coeff,score,latent] = pca(data)

biplot(coeff(:,1:2),'scores',score(:,1:2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 看起来不行，那就更换模拟方法 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x = randn(500,1)*25 + 50
y = randn(500,1)*25 + 50

x(find(x<0),1) = 0
y(find(y<0),1) = 0

save oxygen_glu.mat x y
%% 重新模拟多尺度模型集，不过这次要记录solution


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
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建\Multiscale_Model_ensemble
for i = 1:500
    load(['model-',num2str(i)]);
    model = Model_thermo_;
    r_EX_ = find(contains(model.varNames,'EX_'));
    model.var_lb(r_EX_,:) = 0;

    for ii = 1:length(wait_open_rxn)
        c_ = find(strcmp(model.varNames,['R_',wait_open_rxn{ii,3}]));
        model.var_lb(c_,1) = -100;
    end

    r_oxygen = find(strcmp(model.varNames,'R_rep_EX_o2_e'));
    r_glucose = find(strcmp(model.varNames,'R_rep_EX_glc_e'));
    model.var_ub(r_oxygen,1) = x(i);
    model.var_ub(r_glucose,1) = y(i);
    solution_1 = gurobi_solve_CY(model,'max');
    Res_f{i,1} = model.varNames;
    Res_f{i,2} = solution_1.x;

    clc
    sprintf('正在处理模型%d',i)
end

save 多尺度模型集模拟v2.mat Res_f












