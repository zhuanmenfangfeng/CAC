%% 首先关闭模型中所有的交换反应,并设定相应的培养条件
% cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\1.代谢网络构建\代谢网络draft;
% load Unversal-model_icw773metadd_rxnadd_rxnrep.mat;
% % 确定基本营养物质
% [~,~,basic_substance] = xlsread('biolog_base_composition.csv')
% basic_substance = basic_substance(:,2)
% basic_substance(1,:) = []
% 
% % 由于icw773生物量反应的特殊性，需要额外添加这些物质
% waitaddmets_new = {'cpd00244[e0]';'cpd00104[e0]';'cpd00099[e0]';'cpd00027[e0]'}
% wait_open_rxn = [waitaddmets_new;basic_substance]
% f = @(x)strrep(x,'_e','[e0]')
% wait_open_rxn = cellfun(f,wait_open_rxn,'UniformOutput',false)
% 
% % 确定对应的交换反应
% for i = 1:length(wait_open_rxn)
%     r = find(strcmp(wait_open_rxn{i,1},Universal_Model.mets));
%     c = find(Universal_Model.S(r,:)==-1);
%     for ii = 1:length(c)
%         if length(find(Universal_Model.S(:,c(ii))~=0))==1;
%             wait_open_rxn{i,2} = Universal_Model.rxns{c(ii),1};
%         end
%     end
% end
% 
% for i = 1:size(wait_open_rxn,1)
%     wait_open_rxn{i,3} = strrep(wait_open_rxn{i,2},'(','_');
%     wait_open_rxn{i,3} = strrep(wait_open_rxn{i,3},')','');
% end
% 
% %% 循环内部内容：
% 
% cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建\Multiscale_Model_ensemble
% load model-1.mat
% 
% 
% % 初始化培养条件
% r = find(contains(Model_thermo_.rxns,'_EX_'));
% r_1 = find(contains(Model_thermo_.rxns,'R_'));
% r_2 = intersect(r,r_1);
% EX_rxns = Model_thermo_.rxns(r_2,:);     %% 确定全部交换反应的位置
% 
% for i = 1:length(EX_rxns)
%     r_ = find(strcmp(EX_rxns{i,1},Model_thermo_.varNames));
%     Model_thermo_.var_ub(r_,1) = 0;
% end                                      %% 初始化培养条件
% 
% % 设置相应营养物质的吸收速率
% for i = 1:length(wait_open_rxn)
%     r_pos = find(strcmp(Model_thermo_.varNames,['R_',wait_open_rxn{i,3}]));
% 
% 
% 
% 
% 
% 
% end

%% 发现一个问题，好像可以直接用限制那些本不属于icw773模型的交换反应的通量：
% clear
% load model-1.mat
% r = find(contains(Model_thermo_.rxns,'R_EX_'));
% EX_rxns = Model_thermo_.rxns(r,:);     %% 确定全部交换反应的位置
% 
% for i = 1:length(EX_rxns)
%     r_ = find(strcmp(EX_rxns{i,1},Model_thermo_.varNames));
%     Model_thermo_.var_ub(r_,1) = 0;
% end 

% 检查生物量方程
r = find(Model_thermo_.f==1)
Model_thermo_.varNames(r,:)                     % 目标方程设置正确

% 求解：
gurobi_solve_CY(Model_thermo_,'max')            % 求解正确，接下来就看设置几个目标方程了

% 测试成功，接下来开始运行模型集求解：

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

% new_target = {'F_rep_EX_lys_L_e';'F_rep_EX_pro_L_e';'F_rep_EX_his_L_e'}
% rxn_names_1 = setdiff(rxn_names_1,new_target) 没啥用 本来也不包含交换反应。

% 建立比生长速率模拟的矩阵
growth_mat = []
his_production = []
pro_production = []
lys_production = []

% 建立单反应敲除多约束模型求解算法
for i = 1:500
    cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建\Multiscale_Model_ensemble
    load(['model-',num2str(i)]);
    model = Model_thermo_;

    r = find(contains(model.rxns,'R_EX_'));
    EX_rxns = model.rxns(r,:);     %% 确定全部交换反应的位置
    for ii = 1:length(EX_rxns)
        r_ = find(strcmp(EX_rxns{ii,1},model.varNames));
        model.var_ub(r_,1) = 0;
    end

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
        solution = gurobi_solve_CY(model1,'max');                          % 求解比生长速率
        growth_mat(i,ii) = solution.objval;

        model1.f = zeros(length(model1.f),1);                               % 求解组氨酸比生长速率
        r_ = find(strcmp(model1.varNames,'F_rep_EX_lys_L_e'));
        model1.f(r_,1) = 1;
        solution = gurobi_solve_CY(model1,'max');
        his_production(i,ii) = solution.objval;

        model1.f = zeros(length(model1.f),1);                               % 求解组氨酸比生长速率
        r_ = find(strcmp(model1.varNames,'F_rep_EX_pro_L_e'));
        model1.f(r_,1) = 1;
        solution = gurobi_solve_CY(model1,'max');
        pro_production(i,ii) = solution.objval;

        model1.f = zeros(length(model1.f),1);                               % 求解组氨酸比生长速率
        r_ = find(strcmp(model1.varNames,'F_rep_EX_his_L_e'));
        model1.f(r_,1) = 1;
        solution = gurobi_solve_CY(model1,'max');
        lys_production(i,ii) = solution.objval;

    end
    clc
    sprintf('正在处理模型%d',i)
end
