function [X_Child] = generate_child2(lb, ub,  pop, param)
% this is the DE (differential evolution) algorithm of Rainer Storn, 
% modification based on the implementation of 
% FAC Viana, SURROGATES Toolbox User's Guide, Version 3.0, 2011, available
% at http://sites.google.com/site/felipeacviana/surrogatestoolbox.
% Initial implementation can be found at http://www.icsi.berkeley.edu/~storn/code.html

pop_old = pop.X;
NP = param.popsize;
CR = 0.9;
F = 0.5;
D = size(pop_old, 2);
% index pointer array
ind = randperm(4);
% shuffle locations of vectors
a1  = randperm(NP);
% csvwrite('a1.csv', a1)

lb = repmat(lb, NP, 1);
ub = repmat(ub, NP, 1);
% rotating index array (size NP)
rot = (0:1:NP-1);

% rotate indices by ind(1) positions
rt = rem(rot+ind(1),NP);
% rotate vector locations
a2  = a1(rt+1);

rt = rem(rot+ind(2),NP);
a3  = a2(rt+1);

rt = rem(rot+ind(3),NP);
a4  = a3(rt+1);

rt = rem(rot+ind(4),NP);
a5  = a4(rt+1);


% shuffled population 1
pm1 = pop_old(a1,:);
% shuffled population 2
pm2 = pop_old(a2,:);
% shuffled population 3
pm3 = pop_old(a3,:);
% shuffled population 4
pm4 = pop_old(a4,:);
% shuffled population 5
pm5 = pop_old(a5,:);

% all random numbers < CR are 1, 0 otherwise
mui = rand(NP,D) < CR;
% inverse mask to mui
mpo = mui < 0.5;
    


% differential variation
ui = pm3 + F * (pm1 - pm2);
% crossover
ui = pop_old .* mpo + ui .* mui;
% DE/rand-to-best/1

% correcting violations on the lower bounds of the variables
% these are good to go
maskLB = ui > lb;
% these are good to go
maskUB = ui < ub;
ui     = ui .* maskLB .* maskUB + lb .* (~maskLB) + ub.* (~maskUB);

X_Child = ui;   
end