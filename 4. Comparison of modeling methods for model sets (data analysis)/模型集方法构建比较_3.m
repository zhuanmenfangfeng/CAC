clear
%% 首先登陆两个gap_filling_reaction集合：
load iterative_gapfilling_method.mat
iterative_gapfilling = gap_filled_rnx
load data_alignment_gap_filled_raction.mat
data_alignment = gap_filled_rnx

% 找到他们的代谢途径注释：
[~,~,data_infor] = xlsread('ModelSEED.xls','sheet1');

% 先是iterative_gapfilling：
f1 = @(x)strrep(x,'_c0','')
data_alignment = cellfun(f1,data_alignment,'UniformOutput',false)
for i = 1:size(data_alignment,1)
    if length(find(strcmp(data_alignment{i,1},data_infor(:,11))))>0
        r1 = find(strcmp(data_alignment{i,1},data_infor(:,11)));
        for ii = 1:length(r1)
            if length(data_infor{r1(ii),17})>1
                data_alignment{i,2} = data_infor{r1(ii),17};
                break
            end
        end
    end
end

f2 = @(x)strrep(x,'_REV','')
iterative_gapfilling = cellfun(f2,iterative_gapfilling,'UniformOutput',false)
iterative_gapfilling = unique(iterative_gapfilling)

f3 = @(x)strrep(x,'_c0','');
iterative_gapfilling = cellfun(f3,iterative_gapfilling,'UniformOutput',false);

for i = 1:size(iterative_gapfilling,1)
    if length(find(strcmp(iterative_gapfilling{i,1},data_infor(:,11))))>0
        r1 = find(strcmp(iterative_gapfilling{i,1},data_infor(:,11)));
        for ii = 1:length(r1)
            if length(data_infor{r1(ii),17})>1
                iterative_gapfilling{i,2} = data_infor{r1(ii),17};
                break
            end
        end
    end
end
save pathway_anatation.mat iterative_gapfilling data_alignment

%% 处理数据 可以把反应的频次给列出来（统计每个方法中途径频次出现最高的几项）
f1 = @(x)strrep(x,'KEGG:','')
f2 = @(x)strrep(x,'MetaCyc:','')

load pathway_anatation.mat
Res_data_alignment = {}
for i = 1:size(data_alignment,1)
    if length(data_alignment{i,2})>1
        res = split(data_alignment{i,2},';');
        res = cellfun(f1,res,'UniformOutput',false);
        res = cellfun(f2,res,'UniformOutput',false);
        Res_data_alignment = [Res_data_alignment;res];
    end
end
Res_data_alignment_unique = unique(Res_data_alignment)

for i = 1:length(Res_data_alignment_unique)
    Res_data_alignment_unique{i,2} = sum(strcmp(Res_data_alignment_unique{i,1},Res_data_alignment));
end
Res_data_alignment_unique = sortrows(Res_data_alignment_unique,2,'descend');

% 接下来再来处理迭代填补的部分：

Res_iterative_gapfilling = {}
for i = 1:size(iterative_gapfilling,1)
    if length(iterative_gapfilling{i,2})>1
        res = split(iterative_gapfilling{i,2},';');
        res = cellfun(f1,res,'UniformOutput',false);
        res = cellfun(f2,res,'UniformOutput',false);
        Res_iterative_gapfilling = [Res_iterative_gapfilling;res];
    end
end
Res_iterative_gapfilling_unique = unique(Res_iterative_gapfilling)

for i = 1:length(Res_iterative_gapfilling_unique)
   Res_iterative_gapfilling_unique{i,2} = sum(strcmp(Res_iterative_gapfilling_unique{i,1},Res_iterative_gapfilling));
end
Res_iterative_gapfilling_unique = sortrows(Res_iterative_gapfilling_unique,2,'descend');

% 人工确定要用于制作箱式图的几个途径：
Res_pathway = {' Degradation (Degradation/Utilization/Assimilation)';' Energy-Metabolism (Generation of Precursor Metabolite and Energy)';
    ' AROMATIC-COMPOUNDS-DEGRADATION (Aromatic Compound Degradation)';' Cofactor-Biosynthesis (Cofactor, Prosthetic Group, Electron Carrier, and Vitamin Biosynthesis'
    ;' Carbohydrates-Degradation (Carbohydrate Degradation)';' Gentisate-Degradation (Gentisate Degradation)'}

% 统计词条相应的反应数量以及出现的频率,写个双循环吧，先确定包含哪些反应，随后再确定反应频率。
% 先统计data_alignment
for i = 1:length(Res_pathway)
    res1 = {};
    for ii = 1:size(data_alignment,1)
        if length(data_alignment{ii,2})>1
            if contains(data_alignment{ii,2},Res_pathway{i,1});
                res1 = [res1;data_alignment{ii,1}];
            end
        end
    end
    Res_pathway{i,2} = res1;
end

% 再统计interative_gap_filling
for i = 1:length(Res_pathway)
    res1 = {};
    for ii = 1:size(iterative_gapfilling,1)
        if length(iterative_gapfilling{ii,2})>1
            if contains(iterative_gapfilling{ii,2},Res_pathway{i,1});
                res1 = [res1;iterative_gapfilling{ii,1}];
            end
        end
    end
    Res_pathway{i,3} = res1;
end

% 然后统计每个反应的出现频率：
% 先统计 data_alignment 方法中反应的频率：
clearvars -except Res_pathway

% 先计算data_alignment方法出现的频率
load ensmeble_rxn.mat
Res_data_align = Res_

for i = 1:size(Res_pathway,1)
    if length(Res_pathway{i,2})>1
        rnx_list = Res_pathway{i,2};
        % 统计在模型集中出现的频率
        for ii = 1:size(rnx_list,1)
            b = 0;
            for iii = 1:length(Res_data_align)
                if sum(contains(Res_data_align{iii,1},rnx_list{ii,1}))>0
                    b = b+1;
                else
                    b = b+0;
                end
            end
            rnx_list{ii,2} = b./500;
        end
        Res_pathway{i,4} = rnx_list;
    end
    clc
    sprintf('已经进行了%d',i)
end

% 再计算iterative_gap_filling的方法的出现频率
load iterative_gap_filling.mat
Res_iterative_gapfilling = Res_

for i = 1:size(Res_pathway,1)
    if length(Res_pathway{i,3})>1
        rnx_list = Res_pathway{i,3};
        % 统计在模型集中出现的频率
        for ii = 1:size(rnx_list,1)
            b = 0;
            for iii = 1:length(Res_iterative_gapfilling)
                if sum(contains(Res_iterative_gapfilling{iii,1},rnx_list{ii,1}))>0
                    b = b+1;
                else
                    b = b+0;
                end
            end
            rnx_list{ii,2} = b./500;
        end
        Res_pathway{i,5} = rnx_list;
    end
    clc
    sprintf('已经进行了%d',i)
end

save pathway_anatation_result.mat Res_pathway