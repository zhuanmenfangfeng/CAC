load Unversal-model.mat
[~,~,final_met_relation] = xlsread('Compound细化分解_2.xlsx','最终代谢物对应关系')
%% 将所有反应都通过增加的方式进入模型
init_model = xls2model('初始模型.xlsx')

[~,~,data_1] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对0包含未知代谢物')

% 1对0不包含未知代谢物字符串未匹配
[~,~,data_2] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对0不包含未知代谢物(代谢字符串未匹配)')

% 1对0不包含未知代谢物字符串匹配
[~,~,data_3] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对0不包含未知代谢物(代谢字符串匹配)')

% 1对1
[~,~,data_4] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对1')

% 1对n 相同
[~,~,data_5] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对n相同')

% 1对n 不同
[~,~,data_6] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对n不同(删减全局反应)+手工精炼')

% Exchange_reaction_wait add
[~,~,data_7] = xlsread('model2017反应注释结果-手工精炼.xlsx','Exchange_rxn_waitadd')

r = find(contains(data_1(:,1),'EX_'))
data_1(r,:) = []

r = find(contains(data_2(:,1),'EX_'))
data_2(r,:) = []

r = find(contains(data_3(:,1),'EX_'))
data_3(r,:) = []

r = find(contains(data_4(:,1),'EX_'))
data_4(r,:) = []

r = find(contains(data_5(:,1),'EX_'))
data_5(r,:) = []

r = find(contains(data_6(:,1),'EX_'))
data_6(r,:) = []

for i = 1:length(data_1)
    c = find(strcmp(data_1{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_1{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

for i = 1:length(data_2)
    c = find(strcmp(data_2{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_2{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

for i = 1:length(data_3)
    c = find(strcmp(data_3{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_3{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

for i = 1:length(data_4)
    c = find(strcmp(data_4{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_4{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

for i = 1:length(data_5)
    c = find(strcmp(data_5{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_5{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

for i = 1:length(data_6)
    c = find(strcmp(data_6{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_6{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

for i = 1:length(data_7)
    met_abr = [data_7{i,1},'[e]'];
    r = find(strcmp(init_model.mets,met_abr));
    c = find(init_model.S(r,:)==-1);
    for ii = 1:length(c)
        if length(find(init_model.S(:,c(ii))~=0))==1
            data_7{i,5} = init_model.rxns{c(ii),1};
        end
    end
end
data_7{163,5} = 'EX_biomass';

for i = 1:length(data_7)
    c = find(strcmp(data_7{i,5},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_7{i,5}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

Universal_Model = changeObjective(Universal_Model,'add_CG_biomass cgl ATCC13032')
optimizeCbModel(Universal_Model)                  % 确实把

r = find(~contains(Universal_Model.rxns,'add_'))
Universal_Model = changeRxnBounds(Universal_Model, Universal_Model.rxns(r,:), 0, 'b')
optimizeCbModel(Universal_Model)              %  比生长速率1000是什么鬼？继续测试，大概率是增加过程的问题

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.metNames = {'test_1'};
model.mets = {'test_1'};
model.b = [0];
model.csense = ['E'];
model.metFormulas = {'test'};
model.metCharge = [0];
model.S = [-1];
model.rxns = {'rxn_test1'};
model.lb = [0];
model.ub = [0];
model.c = [0];
model.rules = {'rules'};
model.subsystems = {'subsytem'}
model.grRules = {'grRules'}
model.rxnNames = {'rxn_test1'}
model.rev = [0]
Universal_Model = model

for i = 1:length(data_1)
    c = find(strcmp(data_1{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_1{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

for i = 1:length(data_2)
    c = find(strcmp(data_2{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_2{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

for i = 1:length(data_3)
    c = find(strcmp(data_3{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_3{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

for i = 1:length(data_4)
    c = find(strcmp(data_4{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_4{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

for i = 1:length(data_5)
    c = find(strcmp(data_5{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_5{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

for i = 1:length(data_6)
    c = find(strcmp(data_6{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_6{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

for i = 1:length(data_7)
    met_abr = [data_7{i,1},'[e]'];
    r = find(strcmp(init_model.mets,met_abr));
    c = find(init_model.S(r,:)==-1);
    for ii = 1:length(c)
        if length(find(init_model.S(:,c(ii))~=0))==1
            data_7{i,5} = init_model.rxns{c(ii),1};
        end
    end
end
data_7{163,5} = 'EX_biomass';

for i = 1:length(data_7)
    c = find(strcmp(data_7{i,5},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_7{i,5}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

Universal_Model = changeObjective(Universal_Model,'add_CG_biomass cgl ATCC13032')
optimizeCbModel(Universal_Model)                  % 确实把增加反应删除反应的函数调校正确了！ 再是一边之前的方法：
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
load waitadd.mat
load Unversal-model.mat
[~,~,final_met_relation] = xlsread('Compound细化分解_2.xlsx','最终代谢物对应关系')
% 首先删除1-6中全部的交换反应
init_model = xls2model('初始模型.xlsx')

r = find(contains(data_1(:,1),'EX_'))
data_1(r,:) = []

r = find(contains(data_2(:,1),'EX_'))
data_2(r,:) = []

r = find(contains(data_3(:,1),'EX_'))
data_3(r,:) = []

r = find(contains(data_4(:,1),'EX_'))
data_4(r,:) = []

r = find(contains(data_5(:,1),'EX_'))
data_5(r,:) = []

r = find(contains(data_6(:,1),'EX_'))
data_6(r,:) = []

for i = 1:length(data_1)
    c = find(strcmp(data_1{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_1{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end


%% data2
for i = 1:length(data_2)
    c = find(strcmp(data_2{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);

    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_2{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

%% data3
for i = 1:length(data_3)
    c = find(strcmp(data_3{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);
    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_3{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

%% data4
for i = 1:size(data_4,1)
    c = find(strcmp(data_4{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);
    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_4{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

%% data5
for i = 1:size(data_5,1)
    c = find(strcmp(data_5{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);
    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_5{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

%% data6
for i = 1:size(data_6,1)
    c = find(strcmp(data_6{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);
    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_6{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end


%% data7
% 首先给data7增加一列反应号码,确定对应的反应
for i = 1:length(data_7)
    met_abr = [data_7{i,1},'[e]'];
    r = find(strcmp(init_model.mets,met_abr));
    c = find(init_model.S(r,:)==-1);
    for ii = 1:length(c)
        if length(find(init_model.S(:,c(ii))~=0))==1
            data_7{i,5} = init_model.rxns{c(ii),1};
        end
    end
end
data_7{20,5} = 'EX_biomass';

for i = 1:size(data_7,1)
    c = find(strcmp(data_7{i,5},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);
    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_7{i,5}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end
save Unversal-model_icw773metadd_rxnadd.mat Universal_Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
load waitrepalce.mat
load Unversal-model_icw773metadd_rxnadd.mat
[~,~,final_met_relation] = xlsread('Compound细化分解_2.xlsx','最终代谢物对应关系')
% 首先删除1-6中全部的交换反应
init_model = xls2model('初始模型.xlsx')

r = find(contains(data_2rep(:,1),'EX_'))
data_2rep(r,:) = []

r = find(contains(data_3rep(:,1),'EX_'))
data_3rep(r,:) = []

r = find(contains(data_4rep(:,1),'EX_'))
data_4rep(r,:) = []

r = find(contains(data_5rep(:,1),'EX_'))
data_5rep(r,:) = []

r = find(contains(data_6rep(:,1),'EX_'))
data_6rep(r,:) = []
%% data_2rep
for i = 1:size(data_2rep,1)
    Universal_Model = removeRxns_CY(Universal_Model,[data_2rep{i,4},'_c0']);
    Universal_Model = removeRxns_CY(Universal_Model,[data_2rep{i,4},'_e0']);
    c = find(strcmp(data_2rep{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);
    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        %
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['rep_',data_2rep{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

%% data_3rep
for i = 1:size(data_3rep,1)
    data_3rep{i,4} = strrep(data_3rep{i,4},'&','')
    data_3rep{i,4} = strrep(data_3rep{i,4},' ','')
    Universal_Model = removeRxns_CY(Universal_Model,[data_3rep{i,4},'_c0']);
    Universal_Model = removeRxns_CY(Universal_Model,[data_3rep{i,4},'_e0']);

    c = find(strcmp(data_3rep{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);
    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['rep_',data_3rep{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

%% data_4rep
for i = 1:size(data_4rep,1)
    data_4rep{i,2} = strrep(data_4rep{i,2},'&','');
    data_4rep{i,2} = strrep(data_4rep{i,2},' ','');
    Universal_Model = removeRxns_CY(Universal_Model,[data_4rep{i,2},'_c0']);
    Universal_Model = removeRxns_CY(Universal_Model,[data_4rep{i,2},'_e0']);
    c = find(strcmp(data_4rep{i,1},init_model.rxns));

    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);
    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['rep_',data_4rep{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

%% data_5rep
pat = 'rxn[0-9]{5}'
for i = 1:size(data_5rep,1)
    a = regexp(data_5rep{i,3},pat,'match');
    for ii = 1:length(a)
        Universal_Model = removeRxns_CY(Universal_Model,[a{ii},'_c0']);
        Universal_Model = removeRxns_CY(Universal_Model,[a{ii},'_e0']);
    end
    c = find(strcmp(data_5rep{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);
    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['rep_',data_5rep{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end

%% data_6rep
pat = 'rxn[0-9]{5}'
for i = 1:size(data_6rep,1)
    a = regexp(data_6rep{i,2},pat,'match');
    Universal_Model = removeRxns_CY(Universal_Model,[a{1},'_c0']);
    Universal_Model = removeRxns_CY(Universal_Model,[a{1},'_e0']);

    c = find(strcmp(data_6rep{i,1},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);
    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['rep_',data_6rep{i,1}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end
%% data_7rep

for i = 1:length(data_7rep)
    met_abr = [data_7rep{i,1},'[e]'];
    r = find(strcmp(init_model.mets,met_abr));
    c = find(init_model.S(r,:)==-1);
    for ii = 1:length(c)
        if length(find(init_model.S(:,c(ii))~=0))==1
            data_7rep{i,5} = init_model.rxns{c(ii),1};
        end
    end
end

for i = 1:size(data_7rep,1)
    Universal_Model = removeRxns_CY(Universal_Model,data_7rep{i,3});
    c = find(strcmp(data_7rep{i,5},init_model.rxns));
    met_pos = find(init_model.S(:,c)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,c);
    for ii = 1:size(met_list,1)                                                           % 直接将所有的信息放在met_list中，如果有需要的话，直接从这里提取相应的信息。
        met_list{ii,3} = met_list{ii,1}(end-2:end);
        met_list{ii,3} = strcat(met_list{ii,3}(1:2),'0',']');
        met_list{ii,2} = met_list{ii,1}(1:end-3);
        met_list{ii,4} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,2},final_met_relation(:,1)));        % 
        met_list{ii,5} = final_met_relation{r1,2};
        met_list{ii,6} = [met_list{ii,5},met_list{ii,3}];
        r2 = find(strcmp(met_list{ii,1},init_model.mets));
        met_list{ii,7} = init_model.metNames{r2,1};
        met_list{ii,8} = init_model.metFormulas{r2,1};
        met_list{ii,9} = 0;
    end
    stoichCoeff = met_list(:,4);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['rep_',data_7rep{i,5}],init_model.rxnNames{c,1},met_list,stoichCoeff,init_model.lb(c),init_model.ub(c));
    clc
    sprintf('已经进行了%d',i)
end
save Unversal-model_icw773metadd_rxnadd_rxnrep.mat Universal_Model
%% 求解测试
sum(contains(Universal_Model.rxns,'rep_')) + sum(contains(Universal_Model.rxns,'add_'))  % 的确是1207
% 将除此之前所有反应的通量关闭：
r1 = find(~contains(Universal_Model.rxns,'rep_'))
r2 = find(~contains(Universal_Model.rxns,'add_'))
r = intersect(r1,r2)
Universal_Model = changeRxnBounds(Universal_Model, Universal_Model.rxns(r,:), 0, 'b')
Universal_Model = changeObjective(Universal_Model,'add_CG_biomass cgl ATCC13032')
optimizeCbModel(Universal_Model)              % 比生长速率0.4342 这么多天第一天能正常运行。

%% 不同碳源利用情况测试：


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 检测代谢是否包含%%%%%%%%%%%%%%
load Unversal-model_icw773metadd_rxnadd_rxnrep.mat
[~,~,carbon_source_test] = xlsread('plata_thresholded.csv')
[~,~,basic_substance] = xlsread('biolog_base_composition.csv')
carbon_source_substance = split(carbon_source_test{1,1},'	')
carbon_source_substance(1,:) = []
basic_substance = basic_substance(:,2)
basic_substance(1,:) = []

% 由于icw773生物量反应的特殊性，需要额外添加这些物质
waitaddmets_new = {'cpd00244[e0]';'cpd00104[e0]';'cpd00099[e0]';'cpd00027[e0]'}

% 首先检测基本物质是否存在于模型：
f = @(x)strrep(x,'_e','[e0]')
basic_substance = cellfun(f,basic_substance,'UniformOutput',false)
basic_substance = [basic_substance;waitaddmets_new]

for i = 1:length(basic_substance)
    if sum(strcmp(basic_substance{i,1},Universal_Model.mets))>0
        basic_substance{i,2} = 0;
    else
        basic_substance{i,2} = 1;
    end
end                      % 全都是0，说明都包含这些代谢物：

f = @(x)strrep(x,'_e','[e0]')
carbon_source_substance = cellfun(f,carbon_source_substance,'UniformOutput',false)
for i = 1:length(carbon_source_substance)
    if sum(strcmp(carbon_source_substance{i,1},Universal_Model.mets))>0
        carbon_source_substance{i,2} = 0;
    else
        carbon_source_substance{i,2} = 1;
    end
end     

% 缺乏5种碳源，再看看碳源的情况：
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

% 删除一下模型中没有的代谢物：
r = find(cell2mat(carbon_source_substance(:,2))==1)
carbon_source_substance(r,:) = []

% 检查上述内容对应的代谢物是否存在于全局模型中：
for i = 1:length(carbon_source_substance)
    r = find(strcmp(carbon_source_substance{i,1},Universal_Model.mets));
    c = find(Universal_Model.S(r,:)~=0);
    for ii = 1:length(c)
        length(find(Universal_Model.S(:,c(ii))~=0))==1;
        carbon_source_substance{i,4} = 'Exchange_rxn_exist';
    end
end              % 此时所有交换反应都存在了。

for i = 1:length(basic_substance)
    r = find(strcmp(basic_substance{i,1},Universal_Model.mets));
    c = find(Universal_Model.S(r,:)~=0);
    for ii = 1:length(c)
        length(find(Universal_Model.S(:,c(ii))~=0))==1;
        basic_substance{i,3} = 'Exchange_rxn_exist';
    end
end              % 此时所有交换反应都存在了。

% 首先关闭所有交换反应，首先得确定所有交换反应吧。
r = find(contains(Universal_Model.rxns,'EX_'))
Universal_Model = changeObjective(Universal_Model,'add_CG_biomass cgl ATCC13032')
changeCobraSolver('gurobi')
optimizeCbModel(Universal_Model)
Universal_Model = changeRxnBounds(Universal_Model, Universal_Model.rxns(r,:), 0, 'l')
optimizeCbModel(Universal_Model)

for ii = 1:length(basic_substance)
    r = find(strcmp(basic_substance{i,1},Universal_Model.mets));
    c = find(Universal_Model.S(r,:)~=0);
    for ii = 1:length(c)
        if length(find(Universal_Model.S(:,c(ii))~=0))==1;
            Universal_Model = changeRxnBounds(Universal_Model,Universal_Model.rxns{c(ii),1},-1000,'l');
        end
    end
end

for i = 1:length(carbon_source_substance)
    Universal_Model_test = Universal_Model;
    r = find(strcmp(carbon_source_substance{i,1},Universal_Model_test.mets));
    c = find(Universal_Model_test.S(r,:)~=0);
    for ii = 1:length(c)
        if length(find(Universal_Model_test.S(:,c(ii))~=0))==1;
            Universal_Model_test = changeRxnBounds(Universal_Model_test,Universal_Model_test.rxns{c(ii),1},-1000,'l');
            sol = optimizeCbModel(Universal_Model_test);
            carbon_source_substance{i,5} = sol.f;
        end
    end
end

% 查看在初始模型上的表现：
init_model = xls2model('初始模型.xlsx')
[~,~,m_open] = xlsread('必须打开的交换反应.csv')

exchange_rxns = init_model.rxns(find(contains(init_model.rxns,'EX_')),:)
init_model = changeRxnBounds(init_model,exchange_rxns,0,'l')
optimizeCbModel(init_model)
init_model = changeRxnBounds(init_model,m_open(:,1),-1000,'l')
optimizeCbModel(init_model)









