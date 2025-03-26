function [solution] = gurobi_solve_CY(model,sence)
Model = [];
Model.varnames = model.varNames;
Model.A = sparse(model.A);
Model.rhs = model.rhs;
Model.obj = model.f;
Model.lb = model.var_lb;
Model.ub = model.var_ub;
Model.modelsense = sence;                      
if iscell(model.constraintType)                            % 这里可以调整的参数： ">", '<','='
    Model.sense = cell2mat(model.constraintType)
else
    Model.sense = model.constraintType
end
Model.sense = char(Model.sense);
Model.vtype = cell2mat(model.vartypes);                    % 这里可以调整的参数： "B", 'C','I'    分别对应的是二元变量 连续变量以及整数变量
% params.PoolSearchMode = 1
% params.PoolSolutions = 50
% % params.PoolGap = 0
% solution = gurobi(Model,params);
solution = gurobi(Model);
end

