function [pop] = reduce_pop(pop, popsize)
% Reduction is based on information of solutions that have been completely evaluated

% 
% pop.X = pop.X(1: popsize, :);
% pop.F = pop.F(1: popsize, :);
% if ~isempty(pop.C)
%     pop.C = pop.C(1: popsize, :);
% end
% 
% if ~isempty(pop.A)
%     pop.A = pop.A(1: popsize, :);
% end
% 
% if ~isempty(pop.Mdl)
%     pop.Mdl = pop.Mdl(1: popsize);
% end
% 
% if ~isempty(pop.trgdata)
%     pop.trgdata = pop.trgdata(1: popsize);
% end


fields = fieldnames(pop);
nf = length(fields);
for i = 1:nf
    if iscell(pop.(fields{i})) % deal with all cell array
        if ~isempty(pop.(fields{i}))
            pop.(fields{i}) = pop.(fields{i})(1:popsize);
        end
    else
        if ~isempty(pop.(fields{i}))
            pop.(fields{i}) = pop.(fields{i})(1:popsize, :);
        end
    end
end

return
