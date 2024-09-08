function [a] = Multi_constraint_model_const()
%% 构建多约束模型集
a = 1
clear

load mymodel_modelseed_metCompSymbol_CompartmentData.mat
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


load thermo_data.mat


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

    model_data = getEnzymeCodes(Model_1);
    kcats      = matchKcats(model_data);
    cd D:\博士期间工作\全新框架的提出\NewToolbox;
    % 增加一个CompartmentData的信息
    % 增加一个metSEEDID属性

    model_data.model = prepModelforTFA(model_data.model,ReactionDB,model_data.model.CompartmentData);

    % 热力学数据准备
    if strcmp(model_data.model.thermo_units,'kJ/mol')
        GAS_CONSTANT = 8.314472/1000; % kJ/(K mol)
        error('Not implemented yet!')
    else
        GAS_CONSTANT = 1.9858775/1000; % Kcal/(K mol)
        % value for the bigM in THermo constraints. This will also be the bound value
        bigMtherm = 1e6;
        DGR_lb = -bigMtherm; %kcal/mol
        DGR_ub =  bigMtherm; %kcal/mol
    end
    TEMPERATURE = 298.15; % K
    RT = GAS_CONSTANT*TEMPERATURE;
    bigM = 1e6;

    model_data.model = checkTransport(model_data.model);
    [num_mets_org, num_rxns] = size(model_data.model.S);

    for ii = 1:num_mets_org                               % 替换代谢物缩写中的小、中括号
        newmetname = model_data.model.mets{ii};
        newmetname = strrep(newmetname,'[','_');
        newmetname = strrep(newmetname,']','');
        newmetname = strrep(newmetname,'(','_');
        newmetname = strrep(newmetname,')','');
        model_data.model.mets{ii} = newmetname;
    end

    for ii = 1:num_rxns                                  % 替换反应缩写中的小、中括号
        newrxnname = model_data.model.rxns{ii};
        newrxnname = strrep(newrxnname,'(','_');
        newrxnname = strrep(newrxnname,')','');
        newrxnname = strrep(newrxnname,'[','_');
        newrxnname = strrep(newrxnname,']','');
        model_data.model.rxns{ii} = newrxnname;
    end

    % 可逆性反应拆分
    Fkcat = kcats.forw.kcats;
    Bkcat = kcats.back.kcats;
    kcats = [Fkcat;Bkcat];

    % Predefine ECnumber and uniprots for enzyme model:
    ECnumbers = [model_data.EC_numbers ; model_data.EC_numbers];
    uniprots  = [model_data.uniprots   ; model_data.uniprots];

    % 根据可逆性拆分反应：
    Model_irrev = ModelToIrreversible(model_data.model);
    enzymes = cell(5000,1);
    [m,n]   = size(uniprots);

    data      = load('ProtDatabase.mat');
    swissprot = data.swissprot;
    kegg      = data.kegg;
    total_rxnID = Model_irrev.rxns;
    enzyme = {};

    Model_enzyme = Model_irrev              % 增加酶约束模块
    for ii = 1:length(total_rxnID)
        [Model_enzyme,enzyme] = enzymesupply(Model_enzyme,total_rxnID{ii},enzyme,uniprots(ii,:),kcats(ii,:),swissprot,kegg);
        clc
        sprintf('已经进行了%d',ii)
    end

    % 生成一个met
    % 对蛋白总量的进一步约束（测试成功）
    Model_enzyme = addReaction(Model_enzyme,'prot_pool_exchange',{'prot_pool'},1,false,0,1000);
    fs = 0.015142593;
    Model_enzyme = changeObjective(Model_enzyme,'F_add_CG_biomass cgl ATCC13032');
    P_pos = find(strcmpi(Model_enzyme.rxns,'prot_pool_exchange'));
    Model_enzyme.ub(P_pos) = fs;
    optimizeCbModel(Model_enzyme);                                          % 还真能求解了，牛逼

   % 可以做一个循环，挨个增加一个变
    Model_enzyme_ = Model_enzyme;                                          % 用来合并同工酶反应速率

    for ii = 1:length(Model_enzyme_.rxns)                                  % 这个应该是找所有的同工酶反应，然后把同工酶的反应通量合并在一起了
        if contains(Model_enzyme_.rxns{ii},'No1')
            rxn_str = strrep(Model_enzyme_.rxns{ii},'No1','');
            rxn_str_ = [rxn_str,'No'];
            r = find(contains(Model_enzyme_.rxns,rxn_str_));
            [num_constr,num_vars] = size(Model_enzyme_.S);
            Model_enzyme_.rxns{num_vars+1,1} = rxn_str;
            Model_enzyme_.c(num_vars+1,1)   = 0;
            Model_enzyme_.lb(num_vars+1,1)   = 0;                          % lower bound
            Model_enzyme_.ub(num_vars+1,1)   = 1000;                       % upper bound10
            Model_enzyme_.S(num_constr,num_vars+1) = 0;                    % 增加一列变量以表示新合并的变量。
            Model_enzyme_.mets{num_constr+1,1} = ['flux_constraint_',rxn_str];                 % strcat('DFSEU_',model.mets{i});
            Model_enzyme_.csense  = [Model_enzyme_.csense;'E'];                                % type of constraint: '<', or '>', or '='
            Model_enzyme_.b(num_constr+1) = 0;
            r = r';
            r = [r,num_vars+1];
            Model_enzyme_.S(num_constr+1,r') = [ones(1,length(r)-1),-1];   % 分别在约束中填补对应的系数
        end
        clc
        sprintf('已经继续了%d',ii)
    end

    optimizeCbModel(Model_enzyme_,'max');                                   % 这里好像好了，我再接着往下看看

    % 修改constraint_type
    f = @(x)strrep(x,'E','=');
    Model_enzyme_.csense = arrayfun(f,Model_enzyme_.csense);

    % 修改var_type
    data = num2cell(ones(size(Model_enzyme_.S,2),1));
    f = @(x)num2str(x);
    data = cellfun(f,data);
    f = @(x)strrep(x,'1','C');
    data = arrayfun(f,data);
    data = num2cell(data);

    % 创建变量与约束存放位置
    Model_ = Model_irrev;                                                             % 用来放置热力学约束
    Model_.A = Model_enzyme_.S;
    Model_.constraintType = Model_enzyme_.csense;

    Model_.constraintNames = Model_enzyme_.mets;
    Model_.rhs = Model_enzyme_.b;
    Model_.varNames = Model_enzyme_.rxns;
    Model_.vartypes = data;
    Model_.var_lb = Model_enzyme_.lb;
    Model_.var_ub = Model_enzyme_.ub;
    objective = Model_enzyme_.rxns(find(Model_enzyme_.c));
    Model_.f = zeros(size(Model_.A,2),1);
    Model_.f(find(ismember(Model_.varNames,objective))) = 1;

    for ii = 1:length(Model_.mets)
        metformula = Model_.metFormulas{ii};
        metDeltaGF = Model_.metDeltaGFtr(ii);   %代谢物对应的DeltaGFtr
        metComp = Model_.metCompSymbol{ii};     %代谢物的分室
        Comp_index = find(ismember(Model_.CompartmentData.compSymbolList,metComp));  %处于什么样的分室中；
        metLConc_lb = log(Model_.CompartmentData.compMinConc(Comp_index));  %分室对应的最小代谢物浓度求对数；
        metLConc_ub = log(Model_.CompartmentData.compMaxConc(Comp_index));  %分室对应的最大代谢物浓度求对数；
        Comp_pH = Model_.CompartmentData.pH(Comp_index);   %分室对应的pH。
        if strcmp(metformula,'H2O');
            Model_ = addNewVariable_CY(Model_, strcat('Conc_',Model_.mets{ii}),'C',[0 0]);         %
        elseif strcmp(metformula,'H');
            Model_ = addNewVariable_CY(Model_, strcat('Conc_',Model_.mets{ii}),'C',[log(10^(-Comp_pH)) log(10^(-Comp_pH))]);
        elseif strcmp(Model.metSEEDID{i},'cpd11416')                           %%%修改，改成生物至的化合物号码。
            % we do not create the thermo variables for biomass metabolite
        elseif  (metDeltaGF < 1E6)
            Model_ = addNewVariable_CY(Model_, strcat('Conc_',Model_.mets{ii}),'C',[metLConc_lb metLConc_ub]);
        end
        clc
        sprintf('已经进行了%d',ii)
    end

    total_rxnID = model_data.model.rxns;
    Model_thermo = Model_;
    for ii = 1:length(total_rxnID)
        [Model_thermo] = thermosupply(Model_thermo,total_rxnID{ii},RT,bigM,ReactionDB);
        clc
        sprintf('已经进行了%d',ii)
    end
end

