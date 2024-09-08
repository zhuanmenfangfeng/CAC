cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\1.代谢网络构建\代谢网络draft
changeCobraSolver('gurobi')

[S] = mmread('model_matrix.mtx')
[~,~,reaction_upper_bound] = xlsread('uppbnd.csv')
[~,~,reaction_lower_bound] = xlsread('lowbnd.csv')
[~,~,met_att] = xlsread('met_att.csv')
[~,~,met_comp] = xlsread('met_att.csv')
[~,~,met_id] = xlsread('met_id.csv')
[~,~,met_name] = xlsread('met_name.csv')
[~,~,obj] = xlsread('obj.csv')
[~,~,react_id] = xlsread('react_id.csv')
[~,~,react_name] = xlsread('react_name.csv')
[~,~,react_rev] = xlsread('react_rev.csv')

met_id(1,:) = []
Universal_Model.mets = met_id(:,2)
met_name(1,:) = []
Universal_Model.metNames = met_name(:,2)
met_comp(1,:) = []
% met_comp{8463,2} = 0
% met_comp{8463,3} = '_'
% met_comp{8464,2} = 0
% met_comp{8464,3} = '_'
Universal_Model.metFormulas = met_comp(:,3)
Universal_Model.metCharge = cell2mat(met_comp(:,2))
Universal_Model.genes = {}
Universal_Model.rxnGeneMat = []
Universal_Model.grRules = {}
react_id(1,:) = []
Universal_Model.rxns = react_id(:,2)
react_name(1,:) = []
Universal_Model.rxnNames = react_name(:,2)
Universal_Model.subSystems = {}
Universal_Model.S = S
reaction_lower_bound(1,:) = []
Universal_Model.lb = cell2mat(reaction_lower_bound(:,2))
reaction_upper_bound(1,:) = []
Universal_Model.ub = cell2mat(reaction_upper_bound(:,2))
Universal_Model.b = zeros(length(Universal_Model.mets),1)
obj(1,:) = []
Universal_Model.c = cell2mat(obj(:,2))
react_rev(1,:) = []
Universal_Model.rev = cell2mat(react_rev(:,2))

mat = []
for i = 1:length(Universal_Model.rxnNames)
    if isnan(Universal_Model.rxnNames{i,1})
        Universal_Model.rxnNames{i,1} = 'empty';
    end
    if contains(Universal_Model.rxnNames{i,1},'Exchange');
       mat = [mat;i]
    end
end

for i = 1:length(mat)
    Universal_Model.lb(mat(i)) = -1000;
end

optimizeCbModel(Universal_Model)

Universal_Model.c = cell2mat(obj(:,2))
Universal_Model.c = zeros(length(obj),1)
Universal_Model.c(15001) = 1

optimizeCbModel(Universal_Model)

model2xls(Universal_Model,'Modeltest.xls')
%%
% 失败，还是得把空的地方填补起来：
grRules = num2cell(ones(length(Universal_Model.rxns),1))
Universal_Model.grRules = grRules
Universal_Model.subSystems = grRules

rxnGeneMat = zeros(length(Universal_Model.rxns),1000);
Universal_Model.rxnGeneMat = sparse(rxnGeneMat);
Universal_Model.genes = grRules;
model2xls(Universal_Model,'Success_universe_model.xls')

% 重新导入一个不能求解的模型：

[S] = mmread('model_matrix.mtx')
[~,~,reaction_upper_bound] = xlsread('uppbnd.csv')
[~,~,reaction_lower_bound] = xlsread('lowbnd.csv')
[~,~,met_att] = xlsread('met_att.csv')
[~,~,met_comp] = xlsread('met_att.csv')
[~,~,met_id] = xlsread('met_id.csv')
[~,~,met_name] = xlsread('met_name.csv')
[~,~,obj] = xlsread('obj.csv')
[~,~,react_id] = xlsread('react_id.csv')
[~,~,react_name] = xlsread('react_name.csv')
[~,~,react_rev] = xlsread('react_rev.csv')

met_id(1,:) = []
Universal_Model.mets = met_id(:,2)
met_name(1,:) = []
Universal_Model.metNames = met_name(:,2)
met_comp(1,:) = []
% met_comp{8463,2} = 0
% met_comp{8463,3} = '_'
% met_comp{8464,2} = 0
% met_comp{8464,3} = '_'
Universal_Model.metFormulas = met_comp(:,3)
Universal_Model.metCharge = cell2mat(met_comp(:,2))
Universal_Model.genes = {}
Universal_Model.rxnGeneMat = []
Universal_Model.grRules = {}
react_id(1,:) = []
Universal_Model.rxns = react_id(:,2)
react_name(1,:) = []
Universal_Model.rxnNames = react_name(:,2)
Universal_Model.subSystems = {}
Universal_Model.S = S
reaction_lower_bound(1,:) = []
Universal_Model.lb = cell2mat(reaction_lower_bound(:,2))
reaction_upper_bound(1,:) = []
Universal_Model.ub = cell2mat(reaction_upper_bound(:,2))
Universal_Model.b = zeros(length(Universal_Model.mets),1)
obj(1,:) = []
Universal_Model.c = cell2mat(obj(:,2))
react_rev(1,:) = []
Universal_Model.rev = cell2mat(react_rev(:,2))

mat = []
for i = 1:length(Universal_Model.rxnNames)
    if isnan(Universal_Model.rxnNames{i,1})
        Universal_Model.rxnNames{i,1} = 'empty';
    end
    if contains(Universal_Model.rxnNames{i,1},'Exchange');
       mat = [mat;i]
    end
end

for i = 1:length(mat)
    Universal_Model.lb(mat(i)) = -1000;
end

grRules = num2cell(ones(length(Universal_Model.rxns),1))
Universal_Model.grRules = grRules
Universal_Model.subSystems = grRules

rxnGeneMat = zeros(length(Universal_Model.rxns),1000);
Universal_Model.rxnGeneMat = sparse(rxnGeneMat);
Universal_Model.genes = grRules;
model2xls(Universal_Model,'defeat_universe_model.xls')

%% 测试只保留生物量相关的反应：
Universal_Model = removeRxns(Universal_Model,Universal_Model.rxns(15163:15170))
Universal_Model = changeObjective(Universal_Model,'rxn12985_c0')
optimizeCbModel(Universal_Model)     % 还是有解，说明增加的几个反应除了生物量并不关键

Universal_Model = removeRxns(Universal_Model,Universal_Model.rxns(15164:15165))
optimizeCbModel(Universal_Model)     % 还有解

Universal_Model = removeRxns(Universal_Model,Universal_Model.rxns(15163))
optimizeCbModel(Universal_Model)     % 无解了

%% 生物方程整合
[~,~,data] = xlsread('对应关系.xlsx','Sheet1')
r = find(strcmp(data(:,2),'自定义'))
data_1 = data(r,:)
data_2 = data_1(:,1)
data_3 = unique(data_2)
% 去重
for i = 1:length(data_3)
    data_3{i,2} = strcat('icw_exclusive_met_',num2str(i));
end

%% 向全局模型中增加代谢物
% 确定全局模型
Universal_Model = xls2model('Success_universe_model.xls');
optimizeCbModel(Universal_Model) % 有解
Universal_Model = removeRxns(Universal_Model,Universal_Model.rxns(15163:15170))
optimizeCbModel(Universal_Model)% 有解
Universal_Model = removeRxns(Universal_Model,Universal_Model.rxns(15164))
optimizeCbModel(Universal_Model)% 有解
Universal_Model = removeRxns(Universal_Model,Universal_Model.rxns(15164))
optimizeCbModel(Universal_Model)% 有解
Universal_Model = changeObjective(Universal_Model,'EX_cpd11416_e0')
optimizeCbModel(Universal_Model)% 有解
Universal_Model = removeRxns(Universal_Model,Universal_Model.rxns(15163))
optimizeCbModel(Universal_Model)
save Unversal-model.mat Universal_Model

% 首先向全局反应中增加代谢物
[~,~,meta_infor] = xlsread("初始模型.xlsx","Metabolite List")
meta_infor(1,:) = []

load Unversal-model.mat
for i = 1:length(data_3)
    r = find(strcmp(data_3{i,1},meta_infor(:,1)));
    data_3{i,3} = meta_infor{r,2};
    data_3{i,4} = meta_infor{r,3};
    data_3{i,5} = meta_infor{r,4};
end


for i = 1:size(data_3,1)
    Universal_Model.mets = addMetabolite();

end

Universal_Model.S = full(Universal_Model.S)
for i = 1:size(data_3,1)
    Universal_Model.mets = cat(1,Universal_Model.mets,{[data_3{i,2},'[c0]']});
    Universal_Model.mets = cat(1,Universal_Model.mets,{[data_3{i,2},'[e0]']});
    Universal_Model.metNames = cat(1,Universal_Model.metNames,{[data_3{i,2},'[c0]']});
    Universal_Model.metNames = cat(1,Universal_Model.metNames,{[data_3{i,2},'[e0]']});
    Universal_Model.metFormulas = cat(1,Universal_Model.metFormulas,data_3(i,4));
    Universal_Model.metFormulas = cat(1,Universal_Model.metFormulas,data_3(i,4));
    Universal_Model.metCharge = [Universal_Model.metCharge;data_3{i,5}];
    Universal_Model.metCharge = [Universal_Model.metCharge;data_3{i,5}];
    Universal_Model.b = [Universal_Model.b;0;0];
    Universal_Model.S = [Universal_Model.S;zeros(2,length(Universal_Model.rxns))];
    Universal_Model.csense = [Universal_Model.csense;"E";"E"];
    clc
    sprintf('已经进行了%d',i)
end

%% 向全局模型中增加相应的生物质反应
f = @(x)length(x)
r_1 = find(cellfun(f,data(:,2))==1)
wait_add = data(r_1,1)
% 确定相应的代谢物ID
for i = 1:length(data)
    if ~strcmp(data{i,2},'自定义') && length(data{i,2})>1;
        data{i,3} = data{i,2};
    end
    if strcmp(data{i,2},'自定义')
        r_2 = find(strcmp(data{i,1},data_3(:,1)));
        data{i,3} = data_3{r_2,2};
    end
end

% 确定相应的代谢物对应的分室
model_init = xls2model('初始模型.xlsx')

for i = 1:length(data)
    if length(data{i,2})==1
        rxn_pos = find(strcmp(data{i,1},model_init.rxns));
        met_pos = find(model_init.S(:,rxn_pos) ~= 0);
        met_list = model_init.mets(met_pos,1);
        for ii = 1:size(met_list,1)
            met_list{ii,2} = met_list{ii,1}(end-2:end);
            met_list{ii,1} = met_list{ii,1}(1:end-3);
            met_list{ii,3} = model_init.S(met_pos(ii),rxn_pos)
        end
    else
        r_3 = find(strcmp(data{i,1},met_list(:,1)));
        str = strrep(data{i,3},' ','');
        data{i,4} = [str,met_list{r_3,2}];
        data{i,5} = met_list{r_3,3};
    end
end

% 更换标识符
for i = 1:length(data)
    data{i,4} = strrep(data{i,4},'[c]','[c0]');
    data{i,4} = strrep(data{i,4},'[e]','[e0]');
end

% 增加反应的可逆性、以及反应的上下限等信息
for i = 1:length(data)
       if length(data{i,2})==1
           rxn_pos = find(strcmp(data{i,1},model_init.rxns));
           data{i,3} = 1;
           data{i,4} = Universal_Model.lb(rxn_pos,1);
           data{i,5} = Universal_Model.ub(rxn_pos,1);
       end
end

% 向全局模型增加反应
f = @(x)length(x)
r_4 = find(cellfun(f,data(:,2))==1)
r_4 = [r_4;206]

for i = 1:length(r_4)-1
    Metlist = data(r_4(i)+1:r_4(i+1)-1,4);
    stoichCoeff = data(r_4(i)+1:r_4(i+1)-1,5);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReaction(Universal_Model,data{r_4(i),1},'metaboliteList',Metlist,'stoichCoeffList' ...
        ,stoichCoeff','lowerBound',data{r_4(i),4},'upperBound',data{r_4(i),5});
end

%% 比生长速率求解
% 打开所有交换反应
r_5 = find(contains(Universal_Model.rxns,'EX_'))
exchange_rxn = Universal_Model.rxns(r_5,1)
for i = 1:length(exchange_rxn)
    Universal_Model = changeRxnBounds(Universal_Model,exchange_rxn{i,1},-1000,'l');
end
Universal_Model = changeObjective(Universal_Model,'CG_biomass cgl ATCC13032')
optimizeCbModel(Universal_Model)    % 确实无解，无需把方程补充回模型


save Unversal-model(icw_773_biomass).mat Universal_Model 

%% 将初始模型中的所有反应补充或者替换至全局模型
% 确认初始模型反应条目和全局模型反应条目之间的对应关系
load Unversal-model(icw_773_biomass).mat
already_exist = Universal_Model.rxns(15163:15186)
model_init = xls2model('初始模型.xlsx')
model_waitvalidate = setdiff(model_init.rxns,already_exist)
xlswrite('model2017反应注释结果-手工精炼.xlsx',model_waitvalidate,'等待手工精炼的反应')

% 读取反应间的对应关系
[~,~,infor_1] = xlsread('model2017反应注释结果.xls','Bigg注释UModelseed注释')
[~,~,infor_2] = xlsread('model2017反应注释结果.xls','model2017未被注释')

infor_1(:,3) = []
data_11 = {}
data_1n = {}
data_10 = {}

pat = 'rxn[0-9]{5}'
for i = 1:length(infor_1)
    if length(infor_1{i,2})==1
        data_10 = cat(1,data_10,infor_1(i,:));
    else
        if length(regexp(infor_1{i,2},pat,'match'))==0
            data_10 = cat(1,data_10,infor_1(i,:));
        elseif length(regexp(infor_1{i,2},pat,'match'))==1
            data_11 = cat(1,data_11,infor_1(i,:));
        else
            data_1n = cat(1,data_1n,infor_1(i,:));
        end
    end
end
data_10(:,2) = []
data_10 = [data_10;infor_2]    % 和正好是1207，没毛病

% 写入Excel表格
xlswrite('model2017反应注释结果-手工精炼.xlsx',data_10,'1对0')
xlswrite('model2017反应注释结果-手工精炼.xlsx',data_11,'1对1')
xlswrite('model2017反应注释结果-手工精炼.xlsx',data_1n,'1对n')

%% 进一步缩减1对N
[~,~,data_1n] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对n')
[~,~,Modelseed] = xlsread('ModelSEED.xls','信息精简版')

Res_1 = {}                          % 放反应物一样的
Res_2 = {}                          % 放反应物不同的
pat = 'rxn[0-9]{5}'
for i = 1:length(data_1n)
    info = regexp(data_1n{i,2},pat,'match');
    info = info';
    for ii = 1:length(info)
        r_6 = find(strcmp(info{ii,1},Modelseed(:,4)));
        info{ii,2} = Modelseed{r_6,1};
    end
    metmat = info(:,2);
    if length(unique(metmat)) == 1
        Res_1 = cat(1,Res_1,data_1n(i,:));
    else
        Res_2 = cat(1,Res_2,data_1n(i,:));
    end
end
xlswrite('model2017反应注释结果-手工精炼.xlsx',Res_1,'1对n相同')
xlswrite('model2017反应注释结果-手工精炼.xlsx',Res_2,'1对n不同')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%全局模型中已经增加的反应%%%%%%%%%%%%%%%%%%%%%%%%

load Unversal-model(icw_773_biomass).mat
rxns_exit = Universal_Model.rxns(15163:15186)
[~,~,data_10] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对0')
rxns_unannotate = setdiff(data_10,rxns_exit)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%在全局反应中%%%%%%%%%%%%%%%%%%%%%%%%%%
% 确定代谢物是否覆盖完全
[~,~,data_1] = xlsread('Compound细化分解_2.xlsx','全局模型包含')
[~,~,data_2] = xlsread('Compound细化分解_2.xlsx','全局模型不包含')
[~,~,data_3] = xlsread('Compound细化分解_2.xlsx','未找到对应关系的代谢物')
init_model = xls2model('初始模型.xlsx')
met = init_model.mets
f1 = @(x)strrep(x,'[c]','')
f2 = @(x)strrep(x,'[e]','')
met = cellfun(f1,met,'UniformOutput',false)
met = cellfun(f2,met,'UniformOutput',false)
met = unique(met)
all = unique([data_1(:,1);data_2(:,1);data_3(:,1)])
lost = setdiff(met,all)    % 发现漏掉了no

% 
[~,~,data_1] = xlsread('Compound细化分解_2.xlsx','全局模型包含')
[~,~,data_2] = xlsread('Compound细化分解_2.xlsx','全局模型不包含')
[~,~,data_3] = xlsread('Compound细化分解_2.xlsx','未找到对应关系的代谢物')
f1 = @(x)strrep(x,'[c]','')
f2 = @(x)strrep(x,'[e]','')
met = cellfun(f1,met,'UniformOutput',false)
met = cellfun(f2,met,'UniformOutput',false)
met = unique(met)
all = unique([data_1(:,1);data_2(:,1);data_3(:,1)])
intersect(met,all)   % 完全重叠，可以了
% 给未知代谢物赋值：
for i = 1:length(data_3)
    data_3{i,4} = strcat('icw_exclusive_met_',num2str(i))
end
xlswrite('Compound细化分解_2.xlsx',data_3,'未找到对应关系的代谢物')
%%%%%%%%%%%%%%%%%%%%%%%%% 删除全局模型中不包含 1对n不同 的反应条目%%%%%%%%%%%%
clear
[~,~,data_1n] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对n不同')
[~,~,data_universal] = xlsread('Success_universe_model.xls','Reaction List')

pat = 'rxn[0-9]{5}'
for i = 1:length(data_1n)
    data_1n{i,2} = strrep(data_1n{i,2},' ','');
    match_rxn = regexp(data_1n{i,2},pat,'match');
    for ii= length(match_rxn)
        if sum(strcmp(match_rxn{ii},data_universal(:,1)))==0
           data_1n{i,2} = strrep(data_1n{i,2},match_rxn{ii},'');
        end
    end
end
xlswrite('model2017反应注释结果-手工精炼.xlsx',data_1n,'1对n不同(删减全局反应)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%1对0的反应处理%%%%%%%%%%%%%%%%%%%%%
[~,~,data_1] = xlsread('Compound细化分解_2.xlsx','全局模型包含')
[~,~,data_2] = xlsread('Compound细化分解_2.xlsx','全局模型不包含')
[~,~,data_3] = xlsread('Compound细化分解_2.xlsx','未找到对应关系的代谢物')
[~,~,data_rxn10] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对0')
init_model = xls2model('初始模型.xlsx')

unknown_met = [data_2(:,1);data_3(:,1)]
mets_index = init_model.mets
f1 = @(x)strrep(x,'[c]','')
f2 = @(x)strrep(x,'[e]','')
mets_index = cellfun(f1,mets_index,'UniformOutput',false)
mets_index = cellfun(f2,mets_index,'UniformOutput',false)

metindex = []
for i = 1:length(unknown_met)
    r = find(strcmp(unknown_met{i,1},mets_index));
    metindex = [metindex;r];
end
metindex= unique(metindex)

Res = {}   % 包含未存在于全局模型中的代谢物构,反应条目需要从头添加至模型中。
res = {}   % 包含了部分未知代谢物和已知代谢物，需要人手工确定一下
for i = 1:length(data_rxn10)
    c = find(strcmp(data_rxn10{i,1},init_model.rxns));
    r_1 = find(init_model.S(:,c)~= 0);
    if intersect(r_1,metindex) > 0
        Res = [Res;data_rxn10{i,1}];
    else
        res = [res;data_rxn10{i,1}];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 将反应方程转化为cpd号码数组

[~,~,data_10] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对0不包含未知代谢物')

% 确定对应关系
[~,~,relation_ship] = xlsread('Compound细化分解_2.xlsx','全局模型包含')
relation_ship(:,2) = []
f = @(x)strrep(x,' ','')
relation_ship = cellfun(f,relation_ship,'UniformOutput',false)

% 确定反应方程
[~,~,rxn_formulation] = xlsread('初始模型.xlsx','Reaction List')
rxn_formulation = rxn_formulation(:,[1,3])
f1 = @(x)strrep(x,'[c]','')
f2 = @(x)strrep(x,'[e]','')
rxn_formulation = cellfun(f1,rxn_formulation,'UniformOutput',false)
rxn_formulation = cellfun(f2,rxn_formulation,'UniformOutput',false)
for i = 1:length(data_10)
    r = find(strcmp(data_10{i,1},rxn_formulation(:,1)))
    data_10{i,2} = rxn_formulation{r,2}
end

% 替换为cpd字符串
for i = 1:length(data_10)
    res = split(data_10{i,2},' ');
    r = find(cellfun(@(x)length(x),res)==0);
    res(r,:) = [];
    r = find(strcmp(res,'->')==1);
    res(r,:) = [];
    r = find(strcmp(res,'<-')==1);
    res(r,:) = [];
    r = find(strcmp(res,'<=>')==1);
    res(r,:) = [];
    r = find(strcmp(res,'+')==1);
    res(r,:) = [];
    for ii = 1:length(res)
        r1 = find(strcmp(res{ii,1},relation_ship(:,1)));
        if length(r1) >0
            res{ii,2} = relation_ship{r1,2};
        else
            res{ii,2} = '无';
        end
    end
    r = find(strcmp(res(:,2),'无')==1);
    res(r,:) = [];
    res(:,2) = sort(res(:,2));
    str= '';
    for ii = 1:size(res,1)
        str = strcat(str,';',res{ii,2});
    end
    str(1) = '';
    data_10{i,3} = str;
end

[~,~,universalmodel] = xlsread('ModelSEED.xls','信息精简版');
universalmodel = universalmodel(:,1)

for i = 1:length(data_10)
    r = find(strcmp(data_10{i,3},universalmodel(:,1)));
    if length(r)>=1
        str = '';
        for ii = 1:length(r)
            str = strcat(str,'&',universalmodel{r(ii),4});
        end
        data_10{i,4} = str;
    end
end

data_10 = sortrows(data_10,4)
r = find(cellfun(@(x)length(x),data_10(:,4))==0)
data_10_no = data_10(r,:)
data_10_yes = data_10
data_10_yes(r,:) = []
xlswrite('model2017反应注释结果-手工精炼.xlsx',data_10_no,'1对0不包含未知代谢物(代谢字符串未匹配)')
xlswrite('model2017反应注释结果-手工精炼.xlsx',data_10_yes,'1对0不包含未知代谢物(代谢字符串匹配)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 逐个Sheet梳理一遍，看其中哪个Sheet中反应不存在于全局模型中
[~,~,universa_model_rxns] = xlsread('Success_universe_model.xls','Reaction List')
[~,~,IOdaixiezifuweipipei] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对0不包含未知代谢物(代谢字符串未匹配)')
universal_rxns = universa_model_rxns(:,1)
f1 = @(x)strrep(x,'_c0','')
f2 = @(x)strrep(x,'_e0','')
universal_rxns = cellfun(f1,universal_rxns,'UniformOutput',false)
universal_rxns = cellfun(f2,universal_rxns,'UniformOutput',false)
for i = 1:size(IOdaixiezifuweipipei,1)
    if contains(IOdaixiezifuweipipei{i,4},'rxn')
        if sum(strcmp(IOdaixiezifuweipipei{i,4},universal_rxns))==0
            IOdaixiezifuweipipei{i,5} = '全局模型没有';
        end
    end
end
xlswrite('model2017反应注释结果-手工精炼.xlsx',IOdaixiezifuweipipei,'1对0不包含未知代谢物(代谢字符串未匹配)')

% 逐个 1对0不包含未知代谢物(代谢字符串匹配)
[~,~,IOdaixiezifupipei] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对0不包含未知代谢物(代谢字符串匹配)')
for i = 1:size(IOdaixiezifupipei,1)
    if length(IOdaixiezifupipei{i,4})>1
        if contains(IOdaixiezifupipei{i,4},'rxn')
            obj_ = strrep(IOdaixiezifupipei{i,4},'&','');
            if sum(strcmp(obj_,universal_rxns))==0
                IOdaixiezifupipei{i,5} = '全局模型没有';
            end
        end
    end
end
xlswrite('model2017反应注释结果-手工精炼.xlsx',IOdaixiezifupipei,'1对0不包含未知代谢物(代谢字符串匹配)')

% 逐个 1对1
[~,~,onevone] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对1')
for i = 1:size(onevone,1)
    if length(onevone{i,2})>1
        obj_ = strrep(onevone{i,2},'&','');
        if sum(strcmp(obj_,universal_rxns))==0
            onevone{i,3} = '全局模型没有';
        end
    end
end
xlswrite('model2017反应注释结果-手工精炼.xlsx',onevone,'1对1');

% 逐个 1对n
[~,~,onevmany] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对n相同')
onevmany(:,3) = onevmany(:,2)
pat = 'rxn[0-9]{5}'
for i = 1:size(onevmany,1)
    if length(onevmany{i,2})>1
        obj_ = regexp(onevmany{i,2},pat,'match');
        for ii = 1:length(obj_)
            if sum(strcmp(obj_{ii},universal_rxns))==0
                onevmany{i,3} = strrep(onevmany{i,3},obj_{ii},'');
            end
        end
    end
end
for i = 1:size(onevmany,1)
    if ~contains(onevmany{i,3},'rxn')
       onevmany{i,4} = '全局模型没有';
    end
end
xlswrite('model2017反应注释结果-手工精炼.xlsx',onevmany,'1对n相同');

% 逐个 1对n不同(删减全局反应)+手工精炼
[~,~,onevmany_del_universal_manual_curation] = xlsread('model2017反应注释结果-手工精炼.xlsx','1对n不同(删减全局反应)+手工精炼')
for i = 1:size(onevmany_del_universal_manual_curation,1)
    if length(onevmany_del_universal_manual_curation{i,2})>1
        obj_ = strrep(onevmany_del_universal_manual_curation{i,2},'&','');
        if sum(strcmp(obj_,universal_rxns))==0
            onevmany_del_universal_manual_curation{i,3} = '全局模型没有';
        end
    end
end
xlswrite('model2017反应注释结果-手工精炼.xlsx',onevmany_del_universal_manual_curation,'1对n不同(删减全局反应)+手工精炼');

%% 重新确定一遍全部的交换反应
model_init = xls2model('初始模型.xlsx')
r_ex = find(contains(model_init.rxns,'EX_'))
[r,~] = find(model_init.S(:,r_ex) ~= 0)
exchange_rxns = model_init.rxns(r_ex,:)
exchange_rxns(:,2) = model_init.mets(r,:)
f1 = @(x)strrep(x,'[e]','')
exchange_rxns(:,2) = cellfun(f1,exchange_rxns(:,2),'UniformOutput',false)
f2 = @(x)strrep(x,'[c]','')
exchange_rxns(:,2) = cellfun(f2,exchange_rxns(:,2),'UniformOutput',false)

% 要确定那些代谢物是库中有的，那些代谢物是库中没有的：
[~,~,met_in_universal] = xlsread('Compound细化分解_2.xlsx','全局模型包含');
[~,~,met_in_seed_notin_universal] = xlsread('Compound细化分解_2.xlsx','全局模型不包含');
[~,~,met_notin_seed] = xlsread('Compound细化分解_2.xlsx','未找到对应关系的代谢物');

% 确定每个代谢物对应的缩写：
for i = length(met_in_universal):-1:1
   if sum(strcmp(met_in_universal{i,1},exchange_rxns(:,2)))==0
       met_in_universal(i,:) = [];
   end
end

for i = length(met_in_seed_notin_universal):-1:1
   if sum(strcmp(met_in_seed_notin_universal{i,1},exchange_rxns(:,2)))==0
       met_in_seed_notin_universal(i,:) = [];
   end
end

for i = length(met_notin_seed):-1:1
   if sum(strcmp(met_notin_seed{i,1},exchange_rxns(:,2)))==0
       met_notin_seed(i,:) = [];
   end
end

% 发现代谢物大部分都在全局模型中（161个），有1个在seed中不在全局模型中，一个即不在seed
% 也不在全局模型中（biomass）
% 首先确定这161个是不是都在全局模型中有着对应的交换反应：
met_in_universal(:,3) = met_in_universal(:,2)
for i = 1:length(met_in_universal)
    met_in_universal{i,3} = ['EX_',met_in_universal{i,3},'_e0'];
end

for i = 1:length(met_in_universal)
    if sum(strcmp(met_in_universal{i,3},universa_model_rxns(:,1)))==0
        met_in_universal{i,4} = '不存在于全局模型中';
    end
end               
xlswrite('model2017反应注释结果-手工精炼.xlsx',met_in_universal,'Exchange_rxn_waitadd');

% 现在检查一下反应之和是不是1207，看看反应之和是不是1207
% 1对0包含未知代谢物
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

data_all = [data_1(:,1);data_2(:,1);data_3(:,1);data_4(:,1);data_5(:,1);data_6(:,1);data_7(:,1)]
r = find(contains(data_all,'EX_'))
data_all(r,:) = []                                      % 泪目 1207个,unique一下也是1207个。

% 确定需要添加的反应：
for i = length(data_2):-1:1
    if contains(data_2{i,4},'rxn') && length(data_2{i,5})==1
        data_2(i,:) = [];
    end
end

f = @(x)length(x)
r_1 = find(cellfun(f,data_3(:,4))>1)
r_2 = find(cellfun(f,data_3(:,5))==1)
r = intersect(r_1,r_2)
data_3(r,:) = []

f = @(x)length(x)
r = find(cellfun(f,data_4(:,3))==1)
data_4(r,:) = []

f = @(x)length(x)
r = find(cellfun(f,data_5(:,4))==1)
data_5(r,:) = []

f = @(x)length(x)
r = find(cellfun(f,data_6(:,3))==1)
data_6(r,:) = []

f = @(x)length(x)
r = find(cellfun(f,data_7(:,4))==1)
data_7(r,:) = []

save waitadd.mat data_1 data_2 data_3 data_4 data_5 data_6 data_7

% 确定需要增加至全局模型中的代谢物：
clear
[~,~,met_in_universal] = xlsread('Compound细化分解_2.xlsx','全局模型包含');
[~,~,met_in_seed_notin_universal] = xlsread('Compound细化分解_2.xlsx','全局模型不包含');
[~,~,met_notin_seed] = xlsread('Compound细化分解_2.xlsx','未找到对应关系的代谢物');

modelinit = xls2model('初始模型.xlsx')
f1 = @(x)strrep(x,'[c]','')
f2 = @(x)strrep(x,'[e]','')
met = cellfun(f1,modelinit.mets,'UniformOutput',false)
met = cellfun(f2,met,'UniformOutput',false) 
met = unique(met)                                     % 这也是747个

met_list = [met_in_universal(:,1);met_notin_seed(:,1);met_in_seed_notin_universal(:,1)]
intersect(met,met_list)                   % 交集也是747个

save waitaddmet.mat met_in_universal met_in_seed_notin_universal met_notin_seed
%% 向全局模型中增加代谢物和反应：
clear
load Unversal-model.mat
% 首先增加代谢物(优先处理位于modelseed数据库但不在全局模型中的代谢物)
load waitaddmet.mat
init_model = xls2model('初始模型.xlsx')

% 先将遗漏的代谢物信息补充至全局模型中：
% [~,~,metlist_init] = xlsread('初始模型.xlsx','Metabolite List')
% for i = 1:length(met_in_seed_notin_universal)
%     r = find(strcmp(met_in_seed_notin_universal{i,1},metlist_init(:,1)));
%     met_in_seed_notin_universal{i,4} = metlist_init{r,2};
%     met_in_seed_notin_universal{i,5} = metlist_init{r,3};
%     met_in_seed_notin_universal{i,6} = metlist_init{r,4};
% end
% 
% for i = 1:size(met_in_seed_notin_universal,1)
%     Universal_Model.mets = cat(1,Universal_Model.mets,{[met_in_seed_notin_universal{i,3},'[c0]']});
%     Universal_Model.mets = cat(1,Universal_Model.mets,{[met_in_seed_notin_universal{i,3},'[e0]']});
%     Universal_Model.metNames = cat(1,Universal_Model.metNames,{[met_in_seed_notin_universal{i,4},'[c0]']});
%     Universal_Model.metNames = cat(1,Universal_Model.metNames,{[met_in_seed_notin_universal{i,4},'[e0]']});
%     Universal_Model.metFormulas = cat(1,Universal_Model.metFormulas,met_in_seed_notin_universal{i,5});
%     Universal_Model.metFormulas = cat(1,Universal_Model.metFormulas,met_in_seed_notin_universal{i,5});
%     Universal_Model.metCharge = [Universal_Model.metCharge;met_in_seed_notin_universal{i,6}];
%     Universal_Model.metCharge = [Universal_Model.metCharge;met_in_seed_notin_universal{i,6}];
%     Universal_Model.b = [Universal_Model.b;0;0];
%     Universal_Model.S = [Universal_Model.S;zeros(2,length(Universal_Model.rxns))];
%     Universal_Model.csense = [Universal_Model.csense;"E";"E"];
%     clc
%     sprintf('已经进行了%d',i)
% end
% 
% for i = 1:length(met_notin_seed)
%     r = find(strcmp(met_notin_seed{i,1},metlist_init(:,1)));
%     met_notin_seed{i,5} = metlist_init{r,2};
%     met_notin_seed{i,6} = metlist_init{r,3};
%     met_notin_seed{i,7} = metlist_init{r,4};
% end
% 
% for i = 1:size(met_notin_seed,1)
%     Universal_Model.mets = cat(1,Universal_Model.mets,{[met_notin_seed{i,4},'[c0]']});
%     Universal_Model.mets = cat(1,Universal_Model.mets,{[met_notin_seed{i,4},'[e0]']});
%     Universal_Model.metNames = cat(1,Universal_Model.metNames,{[met_notin_seed{i,5},'[c0]']});
%     Universal_Model.metNames = cat(1,Universal_Model.metNames,{[met_notin_seed{i,5},'[e0]']});
%     Universal_Model.metFormulas = cat(1,Universal_Model.metFormulas,met_notin_seed{i,6});
%     Universal_Model.metFormulas = cat(1,Universal_Model.metFormulas,met_notin_seed{i,6});
%     Universal_Model.metCharge = [Universal_Model.metCharge;met_notin_seed{i,7}];
%     Universal_Model.metCharge = [Universal_Model.metCharge;met_notin_seed{i,7}];
%     Universal_Model.b = [Universal_Model.b;0;0];
%     Universal_Model.S = [Universal_Model.S;zeros(2,length(Universal_Model.rxns))];
%     Universal_Model.csense = [Universal_Model.csense;"E";"E"];
%     clc
%     sprintf('已经进行了%d',i)
% end  
% 
% met_in_seed_notin_universal(:,2) = []
% met_in_universal(:,2) = []
% met_notin_seed(:,2:3) = []
% final_met_relation = [met_in_seed_notin_universal;met_in_universal;met_notin_seed]
% xlswrite('Compound细化分解_2.xlsx',final_met_relation,'最终代谢物对应关系')
% save Unversal-model_icw773metadd.mat Universal_Model

%% 增加反应 
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

%% data5
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

for i = 1:length(data_7)                                                   % 修改
    r = find(strcmp(data_7{i,5},init_model.rxns));                         % 修改
    met_pos = find(init_model.S(:,r)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,r);
    for ii = 1:length(met_list)
        met_list{ii,2} = met_list{ii,1}(end-2:end);
        met_list{ii,2} = strcat(met_list{ii,2}(1:2),'0',']');
        met_list{ii,1} = met_list{ii,1}(1:end-3);
        met_list{ii,3} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,1},final_met_relation(:,1)));        %
        met_list{ii,4} = final_met_relation{r1,2};
        met_list{ii,5} = [met_list{ii,4},met_list{ii,2}];
    end
    Metlist = met_list(:,5);
    stoichCoeff = met_list(:,3);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['add_',data_7{i,5}],init_model.rxnNames{r,1},Metlist,stoichCoeff,init_model.lb(r),init_model.ub(r));
    clc
    sprintf('已经进行了%d',i)
end
save Unversal-model_icw773metadd_rxnadd.mat Universal_Model

%% 替换反应
clear
% 现在检查一下反应之和是不是1207，看看反应之和是不是1207
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

% 确定需要添加的反应：
data_2rep = {}
for i = 1:length(data_2)
    if contains(data_2{i,4},'rxn') && length(data_2{i,5})==1
        data_2rep = cat(1,data_2rep,data_2(i,:));
    end
end
r = find(contains(data_2rep(:,1),'EX_'))
data_2rep(r,:) = []

data_3rep = {}
f = @(x)length(x)
r_1 = find(cellfun(f,data_3(:,4))>1)
r_2 = find(cellfun(f,data_3(:,5))==1)
r = intersect(r_1,r_2)
data_3rep = data_3(r,:)
r = find(contains(data_3rep(:,1),'EX_'))
data_3rep(r,:) = []


data_4rep = {}
f = @(x)length(x)
r = find(cellfun(f,data_4(:,3))==1)
data_4rep = data_4(r,:)
r = find(contains(data_4rep(:,1),'EX_'))
data_4rep(r,:) = []

data_5rep = {}
f = @(x)length(x)
r = find(cellfun(f,data_5(:,4))==1)
data_5rep = data_5(r,:)
r = find(contains(data_5rep(:,1),'EX_'))
data_5rep(r,:) = []

data_6rep = {}
f = @(x)length(x)
r = find(cellfun(f,data_6(:,3))==1)
data_6rep = data_6(r,:)
r = find(contains(data_6rep(:,1),'EX_'))
data_6rep(r,:) = []

data_7rep = {}
f = @(x)length(x)
r = find(cellfun(f,data_7(:,4))==1)
data_7rep = data_7(r,:)
save waitrepalce.mat data_2rep data_3rep data_4rep data_5rep data_6rep data_7rep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% 统一把所有反应都删除掉：

for i = 1:length(data_2rep)                                                % 修改
    Universal_Model = removeRxns_CY(Universal_Model,[data_2rep{i,4},'_c0']);
    Universal_Model = removeRxns_CY(Universal_Model,[data_2rep{i,4},'_e0']);
    r = find(strcmp(data_2rep{i,1},init_model.rxns));                      % 修改
    met_pos = find(init_model.S(:,r)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,r);
    for ii = 1:length(met_list)
        met_list{ii,2} = met_list{ii,1}(end-2:end);
        met_list{ii,2} = strcat(met_list{ii,2}(1:2),'0',']');
        met_list{ii,1} = met_list{ii,1}(1:end-3);
        met_list{ii,3} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,1},final_met_relation(:,1)));        % 
        met_list{ii,4} = final_met_relation{r1,2};
        met_list{ii,5} = [met_list{ii,4},met_list{ii,2}];
    end
    Metlist = met_list(:,5);
    stoichCoeff = met_list(:,3);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['rep_',data_2rep{i,1}],init_model.rxnNames{r,1},Metlist,stoichCoeff,init_model.lb(r),init_model.ub(r));
    clc 
    sprintf('已经进行了%d',i)
end


for i = 1:length(data_3rep)                                                % 修改
    data_3rep{i,4} = strrep(data_3rep{i,4},'&','');                         % 修改
    data_3rep{i,4} = strrep(data_3rep{i,4},' ','');                         % 修改

    Universal_Model = removeRxns_CY(Universal_Model,[data_3rep{i,4},'_c0']);  % 修改
    Universal_Model = removeRxns_CY(Universal_Model,[data_3rep{i,4},'_e0']);  % 修改

    r = find(strcmp(data_4rep{i,1},init_model.rxns));                      % 修改
    met_pos = find(init_model.S(:,r)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,r);
    for ii = 1:length(met_list)
        met_list{ii,2} = met_list{ii,1}(end-2:end);
        met_list{ii,2} = strcat(met_list{ii,2}(1:2),'0',']');
        met_list{ii,1} = met_list{ii,1}(1:end-3);
        met_list{ii,3} = met_coeff(ii,1);                                 
        r1 = find(strcmp(met_list{ii,1},final_met_relation(:,1)));        
        met_list{ii,4} = final_met_relation{r1,2};
        met_list{ii,5} = [met_list{ii,4},met_list{ii,2}];
    end
    Metlist = met_list(:,5);
    stoichCoeff = met_list(:,3);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['rep_',data_3rep{i,1}],init_model.rxnNames{r,1},Metlist,stoichCoeff,init_model.lb(r),init_model.ub(r));
    clc 
    sprintf('已经进行了%d',i)
end


for i = 1:length(data_4rep)                                                % 修改
    data_4rep{i,2} = strrep(data_4rep{i,2},'&','');                         % 修改
    data_4rep{i,2} = strrep(data_4rep{i,2},' ','');                         % 修改
    Universal_Model = removeRxns_CY(Universal_Model,[data_4rep{i,2},'_c0']);  % 修改
    Universal_Model = removeRxns_CY(Universal_Model,[data_4rep{i,2},'_e0']);  % 修改
    r = find(strcmp(data_4rep{i,1},init_model.rxns));                      % 修改
    met_pos = find(init_model.S(:,r)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,r);
    for ii = 1:length(met_list)
        met_list{ii,2} = met_list{ii,1}(end-2:end);
        met_list{ii,2} = strcat(met_list{ii,2}(1:2),'0',']');
        met_list{ii,1} = met_list{ii,1}(1:end-3);
        met_list{ii,3} = met_coeff(ii,1);                                 
        r1 = find(strcmp(met_list{ii,1},final_met_relation(:,1)));        
        met_list{ii,4} = final_met_relation{r1,2};
        met_list{ii,5} = [met_list{ii,4},met_list{ii,2}];
    end
    Metlist = met_list(:,5);
    stoichCoeff = met_list(:,3);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['rep_',data_4rep{i,1}],init_model.rxnNames{r,1},Metlist,stoichCoeff,init_model.lb(r),init_model.ub(r));
    clc 
    sprintf('已经进行了%d',i)
end

pat = 'rxn[0-9]{5}'
for i = 1:length(data_5rep)                                                 % 修改
    a = regexp(data_5rep{i,3},pat,'match');
    data_5rep{i,4} = a{1};
    data_5rep{i,4} = strrep(data_5rep{i,4},'&','');                         % 修改
    data_5rep{i,4} = strrep(data_5rep{i,4},' ','');                         % 修改
    Universal_Model = removeRxns_CY(Universal_Model,[data_5rep{i,4},'_c0']);  % 修改
    Universal_Model = removeRxns_CY(Universal_Model,[data_5rep{i,4},'_e0']);  % 修改
    r = find(strcmp(data_5rep{i,1},init_model.rxns));                      % 修改
    met_pos = find(init_model.S(:,r)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,r);
    for ii = 1:length(met_list)
        met_list{ii,2} = met_list{ii,1}(end-2:end);
        met_list{ii,2} = strcat(met_list{ii,2}(1:2),'0',']');
        met_list{ii,1} = met_list{ii,1}(1:end-3);
        met_list{ii,3} = met_coeff(ii,1);
        r1 = find(strcmp(met_list{ii,1},final_met_relation(:,1)));
        met_list{ii,4} = final_met_relation{r1,2};
        met_list{ii,5} = [met_list{ii,4},met_list{ii,2}];
    end
    Metlist = met_list(:,5);
    stoichCoeff = met_list(:,3);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['rep_',data_5rep{i,1}],init_model.rxnNames{r,1},Metlist,stoichCoeff,init_model.lb(r),init_model.ub(r));
    clc
    sprintf('已经进行了%d',i)
end

pat = 'rxn[0-9]{5}'
for i = 1:size(data_6rep,1)                                                 % 修改
    a = regexp(data_6rep{i,2},pat,'match');
    data_6rep{i,3} = a{1};
    data_6rep{i,3} = strrep(data_6rep{i,3},'&','');                         % 修改
    data_6rep{i,3} = strrep(data_6rep{i,3},' ','');                         % 修改

    Universal_Model = removeRxns_CY(Universal_Model,[data_6rep{i,3},'_c0']);  % 修改
    Universal_Model = removeRxns_CY(Universal_Model,[data_6rep{i,3},'_e0']);  % 修改

    r = find(strcmp(data_6rep{i,1},init_model.rxns));                      % 修改
    met_pos = find(init_model.S(:,r)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,r);
    for ii = 1:length(met_list)
        met_list{ii,2} = met_list{ii,1}(end-2:end);
        met_list{ii,2} = strcat(met_list{ii,2}(1:2),'0',']');
        met_list{ii,1} = met_list{ii,1}(1:end-3);
        met_list{ii,3} = met_coeff(ii,1);                                 
        r1 = find(strcmp(met_list{ii,1},final_met_relation(:,1)));        
        met_list{ii,4} = final_met_relation{r1,2};
        met_list{ii,5} = [met_list{ii,4},met_list{ii,2}];
    end
    Metlist = met_list(:,5);
    stoichCoeff = met_list(:,3);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['rep_',data_6rep{i,1}],init_model.rxnNames{r,1},Metlist,stoichCoeff,init_model.lb(r),init_model.ub(r));
    clc 
    sprintf('已经进行了%d',i)
end

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


for i = 1:length(data_7rep)                                                 % 修改                   % 修改
    Universal_Model = removeRxns_CY(Universal_Model,[data_7rep{i,3},'_c0']);  % 修改
    Universal_Model = removeRxns_CY(Universal_Model,[data_7rep{i,3},'_e0']);  % 修改
    r = find(strcmp(data_7rep{i,5},init_model.rxns));                      % 修改
    met_pos = find(init_model.S(:,r)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,r);
    for ii = 1:length(met_list)
        met_list{ii,2} = met_list{ii,1}(end-2:end);
        met_list{ii,2} = strcat(met_list{ii,2}(1:2),'0',']');
        met_list{ii,1} = met_list{ii,1}(1:end-3);
        met_list{ii,3} = met_coeff(ii,1);                                 
        r1 = find(strcmp(met_list{ii,1},final_met_relation(:,1)));        
        met_list{ii,4} = final_met_relation{r1,2};
        met_list{ii,5} = [met_list{ii,4},met_list{ii,2}];
    end
    Metlist = met_list(:,5);
    stoichCoeff = met_list(:,3);
    stoichCoeff = cell2mat(stoichCoeff);
    Universal_Model = addReactions_CY(Universal_Model,['rep_',data_7rep{i,5}],init_model.rxnNames{r,1},Metlist,stoichCoeff,init_model.lb(r),init_model.ub(r));
    clc 
    sprintf('已经进行了%d',i)
end
save Unversal-model_icw773metadd_rxnadd_rxnrep.mat Universal_Model
%% 求解验证
load Unversal-model_icw773metadd_rxnadd_rxnrep.mat
% 首先检查替换和新增反应之和是不是1207
sum(contains(Universal_Model.rxns,'rep_')) + sum(contains(Universal_Model.rxns,'add_'))  % 的确是1207
% 将除此之前所有反应的通量关闭：
r1 = find(~contains(Universal_Model.rxns,'rep_'))
r2 = find(~contains(Universal_Model.rxns,'add_'))
r = intersect(r1,r2)
Universal_Model = changeRxnBounds(Universal_Model, Universal_Model.rxns(r,:), 0, 'b')
% changeCobraSolver('gurobi')
Universal_Model = changeObjective(Universal_Model,'add_CG_biomass cgl ATCC13032')
optimizeCbModel(Universal_Model)                  % 确实把模型和反应增加进去确实无解，也就是说无解的原因跟全局模型无关，要不是就是代谢物的原因，要不就是我添加模型的方式有问题。

% 检查一下是不是有重复的反应，
rxn_list = Universal_Model.rxns([r1;r2],:)
f1 = @(x)strrep(x,'rep_','')
f2 = @(x)strrep(x,'add_','')
rxn_list = cellfun(f1,rxn_list,'UniformOutput',false)
rxn_list = cellfun(f2,rxn_list,'UniformOutput',false)
rxn_list = unique(rxn_list)                 % 发现没有重复


%% 我就完全用新增反应的方法新建一个模型，查找一下问题：
clear
init_model = xls2model('初始模型.xlsx')
[~,~,final_met_relation] = xlsread('Compound细化分解_2.xlsx','最终代谢物对应关系')

% 首先把代谢物列表拓展一下
model.metNames = {};
model.mets = {};
model.b = [];
model.csense = [];
model.metFormulas = {};
model.metCharge = [];
model.S = [];
model.rxns = {};
model.lb = [];
model.ub = [];
model.c = [];
model.rules = {};
model.subsystems = {}
model.grRules = {}
model.rules = {}
model.rxnNames = {}
model.rev = {}


% 可以先把反应矩阵建立好，最后向矩阵中增加相应的物质信息：

for i = 1:length(final_met_relation)
    metslist = {};
    r = find(strcmp([final_met_relation{i,1},'[c]'],init_model.mets));
    if length(r)==0
        r = find(strcmp([final_met_relation{i,1},'[e]'],init_model.mets));
    end
    metslist{1,1} = final_met_relation{i,2};
    metslist{1,2} = init_model.metNames{r,1};
    metslist{1,3} = init_model.metFormulas{r,1};
    metslist{1,4} = 0;
    model = addMetabolite_CY(model,metslist);
    clc
    sprintf('已经进行了%d',i)
end

model.S = zeros(length(model.mets),length(init_model.rxns))

for i = 1:length(init_model.rxns)                                                   % 修改                     % 修改
    met_pos = find(init_model.S(:,i)~=0);
    met_list = init_model.mets(met_pos,:);
    met_coeff = init_model.S(met_pos,i);
    for ii = 1:size(met_list,1)
        met_list{ii,2} = met_list{ii,1}(end-2:end);
        met_list{ii,2} = strcat(met_list{ii,2}(1:2),'0',']');
        met_list{ii,1} = met_list{ii,1}(1:end-3);
        met_list{ii,3} = met_coeff(ii,1);                                 % 系数
        r1 = find(strcmp(met_list{ii,1},final_met_relation(:,1)));        %
        met_list{ii,4} = final_met_relation{r1,2};
        met_list{ii,5} = [met_list{ii,4},met_list{ii,2}];
    end
    Metlist = met_list(:,5);
    stoichCoeff = met_list(:,3);
    stoichCoeff = cell2mat(stoichCoeff);
    for ii = 1:length(Metlist)
        r = find(strcmp(Metlist{ii,1},model.mets));
        model.S(r,i) = stoichCoeff(ii,1);
    end
    model.rxns = cat(1,model.rxns,{['new_add_',init_model.rxns{i}]});
    model.lb = [model.lb;init_model.lb(i)];
    model.ub = [model.ub;init_model.ub(i)];
    model.c = [model.c;0];
    model.rules = cat(1,model.rules,{'rules'});
    model.subsystems = cat(1,model.subsystems,{'substem'});
    model.grRules = cat(1,model.grRules,{'grRules'});
    model.rxnNames = cat(1,model.rxnNames,{'rxnname'});
    model.rev = [model.rev;0];
    clc
    sprintf('已经进行了%d',i)
end


model = changeObjective(model,'new_add_CG_biomass cgl ATCC13032')
optimizeCbModel(model)                                                     % 我擦 对了0.4342，也就是我这种对应关系是没有问题的，之前的函数和方法存在问题。之前可能还是因为反应上下限的问题，再运行一遍上述内容。





%% Dataadd

Dataadd = [data_1;data_2(:,1);data_3(:,1);data_4(:,1);data_5(:,1);data_6(:,1);data_7(:,5)]
Datarep = [data_2rep(:,1);data_3rep(:,1);data_4rep(:,1);data_5rep(:,1);data_6rep(:,1);data_7rep(:,5)]
data_all = [Dataadd;Datarep]  % 确实好像有重复的内容

f = @(x)sum(strcmp(x,data_all))
r = find(cellfun(f,data_all)>1)
chongfu = unique(data_all(r,:))
icw773_rxn = init_model.rxns
intersect(icw773_rxn,data_all)                        % 反应应该是没问题，重复不了应该也不是这个问题，我现在把所有的反应都增加不替换 进行尝试





































r = find(contains(Universal_Model.rxns,'EX_')==1)
rxnNameList = Universal_Model.rxns(r,:)
Universal_Model = changeRxnBounds(Universal_Model, rxnNameList, 0, 'b')
Universal_Model = changeObjective(Universal_Model,'add_CG_biomass cgl ATCC13032')
changeCobraSolver('gurobi')
optimizeCbModel(Universal_Model)
initCobraToolbox(false)

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


