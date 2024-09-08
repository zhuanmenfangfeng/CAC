%% 就跑一个就行，登录多尺度模型集的模拟结果：
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建\Multiscale_Model_ensemble
load model-1.mat

% 保留热力学与酶约束
model = Model_thermo_;
growthrate = [];
for i = 0:0.1:15
    r_glucose = find(strcmp(model.varNames,'R_rep_EX_glc_e'));
    model.var_ub(r_glucose,1) = i;
    solution = gurobi_solve_CY(model,'max');
    growthrate = [growthrate;solution.objval];
end

glucose_upatake_rate = 0:0.1:15
data= [glucose_upatake_rate',growthrate]
plot(data(:,1),data(:,2))

cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\3.多约束模型集模拟
xlswrite('约束条件的组合与拆分.xlsx',data,'热力学+酶约束');

% 能不能直接在上面的基础上继续分析剩下的三个种模型：
% 酶约束直接把蛋白成分拉到慢就行
% 热力学约束二元变量删了就行：

% lost 酶约束：
model_protein_lost = Model_thermo_
r = find(contains(model_protein_lost.varNames,'prot'))
% 将其上限调整至最大：
model_protein_lost.var_ub(r(end),1) = 10000

% 再模拟比生长速率：
growthrate = [];
for i = 0:0.1:15
    r_glucose = find(strcmp(model_protein_lost.varNames,'R_rep_EX_glc_e'));
    model_protein_lost.var_ub(r_glucose,1) = i;
    solution = gurobi_solve_CY(model_protein_lost,'max');
    growthrate = [growthrate;solution.objval];
end
glucose_upatake_rate = 0:0.1:15
data= [glucose_upatake_rate',growthrate]
plot(data(:,1),data(:,2))
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\3.多约束模型集模拟
xlswrite('约束条件的组合与拆分.xlsx',data,'lost酶约束');


% lost 热力学：
model_thermo_lost = Model_thermo_;
r1 = find(contains(model_thermo_lost.constraintNames,'UF_'));
r2 = find(contains(model_thermo_lost.constraintNames,'UR_'));
r = union(r1,r2);

model_thermo_lost.rhs(r,:) = [];
model_thermo_lost.constraintNames(r,:) = [];
model_thermo_lost.constraintType(r,:) = [];
model_thermo_lost.A(r,:) = [];

% 再模拟比生长速率：
growthrate = [];
for i = 0:0.1:15
    r_glucose = find(strcmp(model_thermo_lost.varNames,'R_rep_EX_glc_e'));
    model_thermo_lost.var_ub(r_glucose,1) = i;
    solution = gurobi_solve_CY(model_thermo_lost,'max');
    growthrate = [growthrate;solution.objval];
end

glucose_upatake_rate = 0:0.1:15
data= [glucose_upatake_rate',growthrate]
plot(data(:,1),data(:,2))
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\3.多约束模型集模拟
xlswrite('约束条件的组合与拆分.xlsx',data,'lost热力学约束');

%% 两个约束都丢失：
% 先丢酶约束
model_protein_thermo_lost = Model_thermo_
r = find(contains(model_protein_thermo_lost.varNames,'pool'))
model_protein_thermo_lost.var_ub(r,1) = 10000

% 再丢热力学约束
r1 = find(contains(model_protein_thermo_lost.constraintNames,'UF_'));
r2 = find(contains(model_protein_thermo_lost.constraintNames,'UR_'));
r = union(r1,r2);

model_protein_thermo_lost.rhs(r,:) = [];
model_protein_thermo_lost.constraintNames(r,:) = [];
model_protein_thermo_lost.constraintType(r,:) = [];
model_protein_thermo_lost.A(r,:) = [];

growthrate = [];
for i = 0:0.1:15
    r_glucose = find(strcmp(model_protein_thermo_lost.varNames,'R_rep_EX_glc_e'));
    model_protein_thermo_lost.var_ub(r_glucose,1) = i;
    solution = gurobi_solve_CY(model_protein_thermo_lost,'max');
    growthrate = [growthrate;solution.objval];
end

glucose_upatake_rate = 0:0.1:15
data= [glucose_upatake_rate',growthrate]
plot(data(:,1),data(:,2))
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\3.多约束模型集模拟
xlswrite('约束条件的组合与拆分.xlsx',data,'lost热力学约束_酶约束');

%% 还是这两个求解一下组氨酸合成途径的吉布斯自由能：
clear
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建\Multiscale_Model_ensemble
load model-1.mat

% 设定好了组氨酸合成途径吉布斯自由能作为目标于此同时设置比生长速率：
model = Model_thermo_;
r_dg = find(contains(model.varNames,'DG_F_rep_HISTD'));
r_his = find(strcmp(model.varNames,'F_rep_HISTD'));
model.f = zeros(length(model.f),1);
model.f(r_dg,1) = 1;
r_growth = find(contains(model.varNames,'F_add_CG_biomass cgl ATCC13032'));

% 把葡萄糖吸收速率设置到10
r_glucose = find(strcmp(model.varNames,'R_rep_EX_glc_e'));
model.var_ub(r_glucose,1) = 10;

% 以组氨酸反应为目标，同时把最大热力学驱动力提出来
mdf = []
his_v = []
for i = 0:0.01:10
    model.var_lb(r_his,1) = i;
    solution = gurobi_solve_CY(model,'min');
    mdf = [mdf;solution.x(r_dg,1)];
end

%% 妈的 热力学驱动力算不出来：

% 直接算算蛋白成本把：

% 多约束模型：
clear
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建\Multiscale_Model_ensemble
load model-1.mat
model = Model_thermo_

% 找到所有的蛋白反应：

r = find(contains(model.varNames,'prot'))
r(end,:) = []

molecular_weight = model.A(:,r)
weight = []
for i = 1:size(molecular_weight,2)
    r = find(molecular_weight(:,i)<0);
    weight = [weight;molecular_weight(r,i)]
end
weight = full(weight)
weight = -weight

r = find(contains(model.varNames,'prot'))
r(end,:) = []

protein_pool = []
for i = 0:0.01:5
    r_glucose = find(strcmp(model.varNames,'R_rep_EX_glc_e'));
    model.var_ub(r_glucose,1) = i;
    solution = gurobi_solve_CY(model,'max');
    protein_pool = [protein_pool;solution.x(r)'*weight];
    clc
    sprintf('已经进行了%d',i)
end

cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\3.多约束模型集模拟
xlswrite('蛋白成本使用情况.xlsx',protein_pool,'多约束模型');


%% 删除热力学再运行一遍：
clear
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\2.多约束模型集的构建\Multiscale_Model_ensemble
load model-1.mat
model = Model_thermo_

r1 = find(contains(model.constraintNames,'UF_'));
r2 = find(contains(model.constraintNames,'UR_'));
r = union(r1,r2);

model.rhs(r,:) = [];
model.constraintNames(r,:) = [];
model.constraintType(r,:) = [];
model.A(r,:) = [];

r = find(contains(model.varNames,'prot'))
r(end,:) = []

molecular_weight = model.A(:,r)
weight = []
for i = 1:size(molecular_weight,2)
    r = find(molecular_weight(:,i)<0);
    weight = [weight;molecular_weight(r,i)]
end
weight = full(weight)
weight = -weight

r = find(contains(model.varNames,'prot'))
r(end,:) = []
r(83) = []

protein_pool = []
weight(83,:) = []
for i = 0:0.01:5
    r_glucose = find(strcmp(model.varNames,'R_rep_EX_glc_e'));
    model.var_ub(r_glucose,1) = i;
    solution = gurobi_solve_CY(model,'max');
    protein_pool = [protein_pool;solution.x(r)'*weight()];
    clc
    sprintf('已经进行了%d',i)
end

cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\3.多约束模型集模拟
xlswrite('蛋白成本使用情况.xlsx',protein_pool,'酶约束模型');

%% 













