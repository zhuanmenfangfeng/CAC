function [] = fastuniversal_construct(Culture_media_path,rxn_score_path,ensemblesize,output_file_path)

[~,~,Culture_media] = xlsread(Culture_media_path,'sheet1')
data_rxn = readtable(rxn_score_path)
Culture_media(1,:) = [];

load Unversal-model_icw773metadd_rxnadd_rxnrep.mat;
Model = Universal_Model;
[~,~,block_reaction] = xlsread('全局模型中的block_reaction.xlsx','Sheet1');
[~,~,repeat_reaction] = xlsread('全局反应中重复的反应.xlsx','Sheet1');
Model = removeRxns(Model,block_reaction);
Model = removeRxns(Model,repeat_reaction);                                 % ok，现在变成了8000多个反应，
wait_remove = {'rxn05746_c0';'rxn12548_c0';'rxn90123_c0'}
Model = removeRxns(Model,wait_remove) 

% 准备五种条件下的矩阵
nutrition_conditions = size(Culture_media,2);
Model.A = sparse(nutrition_conditions.*size(Model.S,1),nutrition_conditions.*size(Model.S,2));
Model.constraintType = {};
Model.constraintNames = {};
Model.rhs = [];

for i = 1:nutrition_conditions
    r_range = [((i-1).*length(Model.mets)+1):i.*length(Model.mets)]';
    c_range = [((i-1).*length(Model.rxns)+1):i.*length(Model.rxns)];
    Model.A(r_range,c_range) = Model.S;                                       
end

b = 1;
% 增加5种培养条件的约束
for i = 1:nutrition_conditions
    for ii=1:size(Model.mets,1)
        Model.constraintType{b,1} = '=';
        Model.constraintNames{b,1} = strcat(num2str(i),'_',Model.mets{ii});  % 代谢物 质量平衡约束
        Model.rhs(b,1) = 0;
        b = b + 1;
    end
end

% 增加5种培养条件的变量
Model.var_lb = [];
Model.var_ub = [];
Model.varNames = {}
Model.vartypes = []
for i = 1:nutrition_conditions
    for ii=1:size(Model.rxns,1)
        Model.var_lb = [Model.var_lb;Model.lb(ii)];
        Model.var_ub = [Model.var_ub;Model.ub(ii)];
        Model.varNames = [Model.varNames;strcat(num2str(i),'_',Model.rxns{ii})];
        Model.vartypes = [Model.vartypes;'C'];
    end
    clc
    sprintf('已经进行了%d',i)
end
Model.f = zeros(size(Model.A,2),1);
Model.vartypes = num2cell(Model.vartypes)

% 增加二元变量，用来限制反应数量
r1 = find(~contains(Model.rxns,'EX_'));
r2 = find(~contains(Model.rxns,'biomass'));
r3 = find(~contains(Model.rxns,'add_'));
r4 = find(~contains(Model.rxns,'rep_'));
r1_2 = intersect(r1,r2);
r1_2_3 = intersect(r1_2,r3);
r1_2_3_4 = intersect(r1_2_3,r4);                                           % 待选取的反应集合
rxn_wait_selected = Model.rxns(r1_2_3_4,:);                                % 试试6000个反应能不能搞出来了

%% 首先为每一种条件都增加一个二元变量用于判断反应是否存在
for i = 1:nutrition_conditions
    for ii = 1:length(rxn_wait_selected)
        Model = addNewVariable_CY(Model, strcat(num2str(i),'_','select_or_not_',rxn_wait_selected{ii}),'B',[0,1]);
        clc
        sprintf('已经进行了%d',ii)
    end
end

% 为每一个反应都增加一个最终变量用于判断反应是否存在
for i = 1:length(rxn_wait_selected)
    Model = addNewVariable_CY(Model, strcat('select_or_not_',rxn_wait_selected{i}),'B',[0,1]);
    clc
    sprintf('已经进行了%d',i)
end

% 对每一种条件都增加相应的约束
for i = 1:nutrition_conditions
    for ii = 1:length(rxn_wait_selected)
        Model.constraintNames = [Model.constraintNames;strcat('Constraint_lb','_',rxn_wait_selected{ii,1})];
        Model.constraintNames = [Model.constraintNames;strcat('Constraint_ub','_',rxn_wait_selected{ii,1})];
        Model.constraintType = [Model.constraintType;{'>'}];
        Model.constraintType = [Model.constraintType;{'<'}];
        pos_rxn = find(strcmp(Model.varNames,strcat(num2str(i),'_',rxn_wait_selected{ii,1})));
        pos_b = find(strcmp(Model.varNames,strcat(num2str(i),'_','select_or_not_',rxn_wait_selected{ii,1})));
        a = sparse(2,size(Model.A,2));
        CLHS.coeffs = [pos_rxn,pos_b];
        a(1,CLHS.coeffs) = [1,-Model.var_lb(pos_rxn,1)];
        a(2,CLHS.coeffs) = [1,-Model.var_ub(pos_rxn,1)];
        Model.A = [Model.A;a];
        Model.rhs = [Model.rhs;0;0];
        clc
        sprintf('已经进行了%d',ii)
    end
end

% 最终增加一个约束用于判断反应是否应该存在
for i = 1:length(rxn_wait_selected)
    pos_b_ii = [];
    for ii = 1:nutrition_conditions
        r_ = find(strcmp(Model.varNames,strcat(num2str(ii),'_select_or_not_',rxn_wait_selected{i})));
        pos_b_ii = [pos_b_ii,r_];
    end
    pos_B = find(strcmp(Model.varNames,strcat('select_or_not_',rxn_wait_selected{i})));
    Model.constraintNames = [Model.constraintNames;strcat('Constraint_final_B',rxn_wait_selected{i,1})];
    Model.constraintType = [Model.constraintType;{'>'}];
    a = sparse(1,size(Model.A,2));
    CLHS.coeffs_1 = [pos_b_ii,pos_B];
    a(1,CLHS.coeffs_1) = [-ones(1,length(pos_b_ii)),1e6];
    Model.A = [Model.A;a];
    Model.rhs = [Model.rhs;0];
    clc
    sprintf('已经进行了%d',i)
end

Model_ = Model

% save big_matrix_example.mat Model_
%% 设置目标函数

data_rxn.normalization = normalize(data_rxn.bitscore,'range')

% r1 = find(~contains(Model.rxns,'EX_'));
% r2 = find(~contains(Model.rxns,'biomass'));
% r3 = find(~contains(Model.rxns,'add_'));
% r4 = find(~contains(Model.rxns,'rep_'));
% r1_2 = intersect(r1,r2);
% r1_2_3 = intersect(r1_2,r3)
% r1_2_3_4 = intersect(r1_2_3,r4);                                           % 待选取的反应集合
% rxn_wait_selected = Model.rxns(r1_2_3_4,:);   

for i = 1:length(rxn_wait_selected)
    rxn_wait_selected{i,2} = 0;
end

f = @(x)strrep(x,'_c0','')
rxn_wait_selected(:,1) = cellfun(f,rxn_wait_selected(:,1),'UniformOutput',false);

% 给反应赋予相应分值：
data_rxn = sortrows(data_rxn,21,'descend');
for i = 1:length(rxn_wait_selected)
    if sum(contains(data_rxn.dbhit,rxn_wait_selected{i,1}))>0
        r_ = find(contains(data_rxn.dbhit,rxn_wait_selected{i,1}));
        rxn_wait_selected{i,2} = data_rxn.normalization(r_(1),1);
        sprintf(rxn_wait_selected{i,1})
    else
        rxn_wait_selected{i,2} = min(data_rxn.normalization);
    end
end
%% 更换不同的培养条件： 
% 初始化培养条件
r_EX_ = find(contains(Model_.varNames,'EX_'));
Model_.var_lb(r_EX_,:) = 0;

% 确定对应的交换反应
Culture_media_Exchange = {}
for i = 1:size(Culture_media,1)
    for ii = 1:size(Culture_media,2)
        r = find(strcmp(Culture_media{i,ii},Universal_Model.mets));
        c = find(Universal_Model.S(r,:)==-1);
        for iii = 1:length(c)
            if length(find(Universal_Model.S(:,c(iii))~=0))==1;
                Culture_media_Exchange{i,ii} = Universal_Model.rxns{c(iii),1};
            end
        end
    end
end

% 设置基本物质的吸收速率
for i = 1:nutrition_conditions
    for ii = 1:size(Culture_media_Exchange,1)
        if length(Culture_media_Exchange{ii,i})>1
            c_ = find(strcmp(Model_.varNames,strcat(num2str(i),'_',Culture_media_Exchange{ii,i})));
            Model_.var_lb(c_,1) = -100;
        end
    end
end



%% 求解
for i = 1:nutrition_conditions
    c_ = find(strcmp(strcat(num2str(i),'_','add_CG_biomass cgl ATCC13032'),Model_.varNames));
    Model_.var_lb(c_,1) = 5;
end
Model_.f = zeros(size(Model_.A,2),1);

f = @(x)strrep(x,'_c0','')
Model_.varNames = cellfun(f,Model_.varNames,'UniformOutput',false);

for i = 1:length(rxn_wait_selected)
    str_1 = ['select_or_not_',rxn_wait_selected{i}];
    r__ = find(strcmp(Model_.varNames,str_1));
    Model_.f(r__,1) = 1./(1+rxn_wait_selected{i,2});    % 两者之前的计算量并没有过多的差异啊，尴尬。
end


% 求解500个试试
solution = gurobi_solve_CY(Model_,'min',ensemblesize);      
r_ = find(solution.x~=0);
Res = [Model_.varNames(r_,1),num2cell(solution.x(r_,1))];

Res_variable_reaction = {};
for i = 1:length(solution.pool)
    r_name = find(solution.pool(i).xn==1);
    for ii = 1:length(r_name)
        Res_variable_reaction{ii,i} = Model_.varNames{r_name(ii),1};
    end
end
xlswrite(output_file_path,Res_variable_reaction,'Sheet1');

end

