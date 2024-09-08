%% 先给每个反应都打上相应的分数
load recorde_rxn_score.mat
data_rxn.normalization = normalize(data_rxn.bitscore,'range');
data_rxn = sortrows(data_rxn,21,'descend');

% 然后挨个找模型，并给相应的反应打上分数，最后在一个循环中给模型打上相应的分数。
% 首先是 data_alignment方法：
load ensmeble_rxn.mat
load data_alignment_gap_filled_raction.mat
f = @(x)strrep(x,'_c0','');

for i = 1:length(Res_)
    data_rxn_ = Res_{i,1};
    % rep和add不参与打分：
    r1 = find(contains(data_rxn_,'add_'));
    r2 = find(contains(data_rxn_,'rep_'));
    r3 = find(contains(data_rxn_,'EX_'));
    r4 = [r1;r2;r3];
    r4 = unique(r4);
    data_rxn_ = setdiff(data_rxn_,data_rxn_(r4,:));
    data_rxn_ = cellfun(f,data_rxn_,'UniformOutput',false);

    for ii = 1:size(data_rxn_,1)
        if sum(contains(data_rxn.dbhit,data_rxn_{ii,1}))>0
            r_ = find(contains(data_rxn.dbhit,data_rxn_{ii,1}));
            data_rxn_{ii,2} = data_rxn.normalization(r_(1),1);
            sprintf(data_rxn_{ii,1});
        else
            data_rxn_{ii,2} = min(data_rxn.normalization);
        end
    end
    s = 0;
    for ii = 1:size(data_rxn_,1)
        s = s + (1./(1 + data_rxn_{ii,2}));
    end
    Res_{i,2} = s;
    clc
    sprintf('已经进行了%d',i)
end

% 开始计算iterative-gap-filling-approch
clearvars -except Res_ data_rxn
data_alignment_score = Res_
load iterative_gap_filling.mat
data_iterative_score = Res_

f = @(x)strrep(x,'_c0','');
f1 = @(x)strrep(x,'_REV','');

for i = 1:length(data_iterative_score)
    data_rxn_ = data_iterative_score{i,1};
    data_rxn_ = cellfun(f,data_rxn_,'UniformOutput',false);
    data_rxn_ = cellfun(f1,data_rxn_,'UniformOutput',false);
    data_rxn_ = unique(data_rxn_);

    for ii = 1:size(data_rxn_,1)
        if sum(contains(data_rxn.dbhit,data_rxn_{ii,1}))>0
            r_ = find(contains(data_rxn.dbhit,data_rxn_{ii,1}));
            data_rxn_{ii,2} = data_rxn.normalization(r_(1),1);
            sprintf(data_rxn_{ii,1});
        else
            data_rxn_{ii,2} = min(data_rxn.normalization);
        end
    end
    s = 0;
    for ii = 1:size(data_rxn_,1)
        s = s + (1./(1 + data_rxn_{ii,2}));
    end
    data_iterative_score{i,2} = s;
    clc
    sprintf('已经进行了%d',i)
end

data_alignment_score(:,2) = num2cell(-cell2mat(data_alignment_score(:,2)))
data_iterative_score(:,2) = num2cell(-cell2mat(data_iterative_score(:,2)))

%% 将结果导出：
xlswrite('模型得分结果.xlsx',data_alignment_score,'data_alignment');
xlswrite('模型得分结果.xlsx',data_iterative_score,'data_iterative');
















