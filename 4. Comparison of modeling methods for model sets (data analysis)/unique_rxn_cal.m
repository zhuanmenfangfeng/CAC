function [num] = unique_rxn_cal(a,Res,solution.pool,r1,rxn_wait_selected,must_need_rxns)
Res_ = {}
for i = 1:length(a)
    Res_1 = Res;
    Res_1(:,3) = num2cell(solution.pool(a(i)).xn);
    r2 = find(cell2mat(Res_1(:,3))==1);
    r = intersect(r1,r2);
    must_need_rxns = Res_1(r,1);
    must_need_rxns = cellfun(f2,must_need_rxns,'UniformOutput',false);
    wait_remove = setdiff(rxn_wait_selected,must_need_rxns);
    Model_1 = removeRxns(Model,wait_remove);
    rxns_list = Model_1.rxns;
    Res_{i,1} = a(i);
    Res_{i,2} = rxns_list;
end
% 统计共同反应：




end

