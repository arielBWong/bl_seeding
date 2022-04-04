function [pop] = reduce_pop(pop, popsize)
% Reduction is based on information of solutions that have been completely evaluated

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
