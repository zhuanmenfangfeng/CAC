clear
cd D:\BaiduSyncdisk\20231009论文写作\数据梳理\5.经典代谢网络模型集和多尺度模型集模拟的结果的区别
load 经典代谢网络模型集模拟v2.mat
load 多尺度模型集模拟v2.mat
% 找到开头是 F_，R_ 的并且结尾不带有No的反应条目汇总一下：

% 多尺度模型这里要合并一下，先统计共涉及哪些反应：
rxns = {}
for i = 1:500
    rxns = [rxns;Res_f_classical{i,1}];
    rxns = unique(rxns);
end

% 检查一下 F_ 和 R_ 是不是都分别存在：
rxns_m = {}
for i = 1:500
    rxns_m = [rxns_m;Res_f{i,1}];
    rxns_m = unique(rxns_m);
end

for i = 1:length(rxns)
    if sum(strcmp(['F_',rxns{i,1}],rxns_m))>0 && sum(strcmp(['R_',rxns{i,1}],rxns_m))>0
        rxns{i,2} = 1;
    end
end
% 检查一下有没有不是0的

sum(cell2mat(rxns(:,2))==0) % 没有

% 那就直接按照rxns的这种方式来取通量就行，没有就赋值为0，先取经典代谢网络模型集：
mat_classical = []
for i = 1:500
    obj_table = [Res_f_classical{i,1},num2cell(Res_f_classical{i,2})];
    for ii = 1:length(rxns)
        if sum(strcmp(rxns{ii,1},obj_table(:,1)))>0
            r = find(strcmp(rxns{ii,1},obj_table(:,1)));
            mat_classical(i,ii) = obj_table{r,2};
        else
            mat_classical(i,ii) = 0;
        end
    end
    clc
    sprintf('已经进行了%d',i)
end

% 再找一下多尺度模型集：
mat_multiscale = []
for i = 1:500
    obj_table = [Res_f{i,1},num2cell(Res_f{i,2})];
    for ii = 1:length(rxns)
        if sum(strcmp(['F_',rxns{ii,1}],obj_table(:,1)))>0 & sum(strcmp(['R_',rxns{ii,1}],obj_table(:,1)))>0
            r1 = find(strcmp(['F_',rxns{ii,1}],obj_table(:,1)));
            r2 = find(strcmp(['R_',rxns{ii,1}],obj_table(:,1)));
            mat_multiscale(i,ii) = obj_table{r1,2}-obj_table{r2,2};
        else
            mat_multiscale(i,ii) = 0;
        end
    end
    clc
    sprintf('已经进行了%d',i)
end

data = [mat_classical;mat_multiscale];

% 删掉全是0的列
[~,c] = find(sum(data)==0);
data(:,c) = [];

[coeff,score,latent] = pca(data);
biplot(coeff(:,1:2),'scores',score(:,1:2));

biplot(coeff(:,1:2),'Scores',score(1:500,1:2),'Color','r','Marker','o');
hold on
biplot(coeff(:,1:2),'Scores',score(501:1000,1:2),'Color','b','Marker','o');









