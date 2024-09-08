%% 生成经典代谢网络模型集
clear
cd D:\博士期间工作\全新框架的提出\NewToolbox
load mymodel_modelseed_metCompSymbol_CompartmentData.mat
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建

load 500_ensemble.mat
load big_matrix_example.mat
load Model_GRP_prepare.mat
Model.CompartmentData = mymodel.CompartmentData
Model.metSEEDID = Model.mets
f1 = @(x)strrep(x,'[c0]','')
f2 = @(x)strrep(x,'[e0]','')
Model.metSEEDID = cellfun(f1,Model.metSEEDID,'UniformOutput',false)
Model.metSEEDID = cellfun(f2,Model.metSEEDID,'UniformOutput',false)
% 不在全局模型中代谢物物直接用NA表示吧
for i = 1:length(Model.metSEEDID)
    if length(Model.metSEEDID{i,1})>8
       Model.metSEEDID{i,1} = 'NA'
    end
end

f1 = @(x)x(1:5);
f2 = @(x)strrep(x,'select_or_not_','');

r1 = find(~contains(Model.rxns,'EX_'));
r2 = find(~contains(Model.rxns,'biomass'));
r3 = find(~contains(Model.rxns,'add_'));
r4 = find(~contains(Model.rxns,'rep_'));
r1_2 = intersect(r1,r2);
r1_2_3 = intersect(r1_2,r3);
r1_2_3_4 = intersect(r1_2_3,r4);                                           % 待选取的反应集合
rxn_wait_selected = Model.rxns(r1_2_3_4,:);                                % 试试6000个反应能不能搞出来了，

Res = Model_.varNames;
Res(:,2) = cellfun(f1,Res(:,1),'UniformOutput',false);
r1 = find(strcmp(Res(:,2),'selec'))                                        % 所有决策变量的位置


Model.metCompSymbol = cellfun(@(x)x(end-3:end),Model.mets,'UniformOutput',false);
Model.metCompSymbol = cellfun(@(x)strrep(x,'[',''),Model.metCompSymbol,'UniformOutput',false);
Model.metCompSymbol = cellfun(@(x)strrep(x,']',''),Model.metCompSymbol,'UniformOutput',false);
Model.metCompSymbol = cellfun(@(x)strrep(x,'0',''),Model.metCompSymbol,'UniformOutput',false);
Model.metCompSymbol = cellfun(@(x)strrep(x,'0',''),Model.metCompSymbol,'UniformOutput',false);

% 向模型中增加一个metComps，胞内是1，胞内是2。
Model.metComps = Model.metCompSymbol;
for i = 1:length(Model.metComps)
    if strcmp(Model.metComps{i,1},'c')
        Model.metComps{i,1} = 1;
    else
        Model.metComps{i,1} = 2;
    end
end
Model.metComps = cell2mat(Model.metComps);

for i = 1:length(solution.pool)

    clearvars -except Model r1 ReactionDB Res rxn_wait_selected i solution f2
    Res_1 = Res;
    Res_1(:,3) = num2cell(solution.pool(i).xn);

    r2 = find(cell2mat(Res_1(:,3))==1);
    r = intersect(r1,r2);                                                  % 所有决策变量不为0的位置

    must_need_rxns = Res_1(r,1);
    must_need_rxns = cellfun(f2,must_need_rxns,'UniformOutput',false);
    wait_remove = setdiff(rxn_wait_selected,must_need_rxns);
    Model_1 = removeRxns(Model,wait_remove);                               % 这就形成了配置好的模型，随后在此基础上进一步组装多尺度模型。

    % 组装多尺度模型，又是一道坎，尽力而为吧，关键首先要考虑适配性的问题。
    % 由于removeRxns在删除反应时不会删除对应的grRules中的内容，根据rules的属性生成相应的grRules
    pat = 'x\([0-9]{1,3}\)'
    Model_1.grRules = {};
    for ii = 1:length(Model_1.rules)
        if length(Model_1.rules{ii,1})>0
            gene_set = regexp(Model_1.rules{ii,1},pat,'match');
            gene_set = gene_set';
            for iii = 1:size(gene_set,1)
                gene_set{iii,2} = Model_1.genes{str2num(gene_set{iii,1}(3:end-1)),1};
            end
            gpr_str = Model_1.rules{ii,1};
            for iii = 1:size(gene_set,1)
                gpr_str = strrep(gpr_str,gene_set{iii,1},gene_set{iii,2});
            end
            Model_1.grRules{ii,1} = gpr_str;
            clc
            sprintf('%d',ii)
        else
            Model_1.grRules{ii,1} = '';
        end

        % 再替换一下字符串，将其中的内容替换成正常的布尔逻辑：
        Model_1.grRules{ii,1} = strrep(Model_1.grRules{ii,1},' | ',' or ');
        Model_1.grRules{ii,1} = strrep(Model_1.grRules{ii,1},' & ',' and ');
    end

    path = ['D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建\Model_ensemble\model-',num2str(i),'.mat'];
    file = 'Model_1';
    save(path,file)
    clc
    sprintf("已经进行了%d",i)
end

%% 经典代谢网络模型集的模拟
clear
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\1.代谢网络构建\代谢网络draft
load Unversal-model_icw773metadd_rxnadd_rxnrep.mat;

[~,~,basic_substance] = xlsread('biolog_base_composition.csv')
basic_substance = basic_substance(:,2)
basic_substance(1,:) = []

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

cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建\Model_ensemble
Res = {};
for i = 1:10
    load(['model-',num2str(i)]);

    Model = Model_1
    % 得想想怎么设置交换反应了
    r_EX_ = find(contains(Model.rxns,'EX_'));
    Model.lb(r_EX_,:) = 0;
    for ii = 1:length(wait_open_rxn)
        c_ = find(strcmp(Model.rxns,wait_open_rxn{ii,2}));
        Model.lb(c_,1) = -100;
    end

    r_oxygen = find(strcmp(Model.rxns,'rep_EX_o2(e)'));
    r_glucose = find(strcmp(Model.rxns,'rep_EX_glc(e)'));

    b = 1
    for ii = 0:3:12
        for iii = 0:3:12
            Model = changeRxnBounds(Model, 'rep_EX_o2(e)', -ii, 'b');
            Model = changeRxnBounds(Model, 'rep_EX_glc(e)', -iii, 'b');
            solution = optimizeCbModel(Model);
            Res{i,b} = solution.f;
            b = b + 1
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 前面的有点问题，现在从solution开始构建模型集：

% 验证一下对不同培养条件的利用能力：
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
Res_f = {};
for i = 1:length(solution.pool)
    Res_1 = Res;
    Res_1(:,3) = num2cell(solution.pool(i).xn);
    r2 = find(cell2mat(Res_1(:,3))==1);
    r = intersect(r1,r2);

    must_need_rxns = Res_1(r,1);
    must_need_rxns = cellfun(f2,must_need_rxns,'UniformOutput',false);
    wait_remove = setdiff(rxn_wait_selected,must_need_rxns);
    Model_1 = removeRxns(Model,wait_remove);
    r_oxygen = find(strcmp(Model_1.rxns,'rep_EX_o2(e)'));
    r_glucose = find(strcmp(Model_1.rxns,'rep_EX_glc(e)'));


    b = 1;
    for ii = 0:3:12
        for iii = 0:3:12
            Model_1 = changeRxnBounds(Model_1, 'rep_EX_o2(e)', -ii, 'l');
            Model_1 = changeRxnBounds(Model_1, 'rep_EX_glc(e)', -iii, 'l');
            solution_1 = optimizeCbModel(Model_1);
            Res_f{i,b} = solution_1.f;
            b = b + 1;
        end
    end
    clc
    sprintf('%d',i)
end                                                        % 这种求解方法获得的解集，90%往上的条件都能够满足。

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 再换个模拟方法，对交换反应空间进行随机取样%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

Res_f_classical = {};
for i = 1:length(solution.pool)
    Res_1 = Res;
    Res_1(:,3) = num2cell(solution.pool(i).xn);
    r2 = find(cell2mat(Res_1(:,3))==1);
    r = intersect(r1,r2);

    must_need_rxns = Res_1(r,1);
    must_need_rxns = cellfun(f2,must_need_rxns,'UniformOutput',false);
    wait_remove = setdiff(rxn_wait_selected,must_need_rxns);
    Model_1 = removeRxns(Model,wait_remove);
    r_oxygen = find(strcmp(Model_1.rxns,'rep_EX_o2(e)'));
    r_glucose = find(strcmp(Model_1.rxns,'rep_EX_glc(e)'));

    Model_1 = changeRxnBounds(Model_1, 'rep_EX_o2(e)', -x(i), 'l');
    Model_1 = changeRxnBounds(Model_1, 'rep_EX_glc(e)', -y(i), 'l');
    solution_1 = optimizeCbModel(Model_1);
    Res_f_classical{i,1} = Model_1.rxns
    Res_f_classical{i,2} = solution_1.x;

    clc
    sprintf('%d',i)
end

save 经典代谢网络模型集模拟v2.mat Res_f_classical












