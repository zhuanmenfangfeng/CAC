function [model] = removeRxns_CY(model,rxns)
inters_rxns = intersect(model.rxns,rxns)
if length(inters_rxns)>0;
    c = []
    for i = 1:length(inters_rxns)
        r_ = find(strcmp(inters_rxns{i,1},model.rxns));
        c = [c;r_];
    end
    model.rxns(c,:) = [];
    model.S(:,c) = [];
    model.lb(c,:) = [];
    model.ub(c,:) = [];
    model.c(c,:) = [];
    model.rxnNames(c,:) = [];
    model.rev(c,:) = [];
    model.rules(c,:) = [];
end
end