function [] = startingPoint_check(seed_xl, prob, xu, match_xl, seed_xu, seed_id, correctionFunc_decision, varargin)
% this function is to find whether local search starting point 
% is in other basin in terms of distance to basic center
% if in other basin show plot

p = inputParser;
addRequired(p,  'seed_xl');
addRequired(p,  'prob');
addRequired(p,  'xu');
addRequired(p,  'match_xl');
addRequired(p,  'seed_xu');
addRequired(p,  'seed_id');
addRequired(p,  'correctionFunc_decision');
addParameter(p, 'delta_mdl', []);
addParameter(p, 'close_mdl', []);
addParameter(p, 'corr_xl', []);
addParameter(p, 'cokrg_mdl', []);
addParameter(p, 'cokrg_xl', []);
addParameter(p, 'trainxl', []);
addParameter(p, 'trainfl', []);
addParameter(p, 'gpr_mdl', []);

% 'co_cheapx',x_trg{1,:}, 'co_cheapy',y_trg{1,:} , 'co_expensivex',x_trg{2,:}  ,'co_expensivey', y_trg{2,:} );
addParameter(p, 'co_cheapx', []);
addParameter(p, 'co_cheapy', []);
addParameter(p, 'co_expensivex', []);
addParameter(p, 'co_expensivey', []);

addParameter(p, 'top_xl', []);

parse(p, seed_xl, prob, xu, match_xl, seed_xu, seed_id, correctionFunc_decision, varargin{:});


seed_xl = p.Results.seed_xl;
prob = p.Results.prob;
xu = p.Results.xu;
match_xl = p.Results.match_xl;
seed_xu = p.Results.seed_xu;
seed_id = p.Results.seed_id;
correctionFunc_decision = p.Results.correctionFunc_decision;
delta_mdl = p.Results.delta_mdl;
close_mdl = p.Results.close_mdl;
corr_xl = p.Results.corr_xl;
cokrg_mdl = p.Results.cokrg_mdl;
cokrg_xl = p.Results.cokrg_xl;
trainxl = p.Results.trainxl;
trainfl = p.Results.trainfl;
co_cheapx = p.Results.co_cheapx;
co_cheapy = p.Results.co_cheapy;
co_expensivex = p.Results.co_expensivex;
co_expensivey = p.Results.co_expensivey;
gpr_mdl = p.Results.gpr_mdl; 
top_xl = p.Results.top_xl;
%--------------

other_centers = mp_othercenters(prob, xu);
xl_prime = prob.get_xlprime(xu);
d_2prime = pdist2(seed_xl, xl_prime);
nc = size(other_centers, 1);

% for i = 1:nc
    %d_2other = pdist2(seed_xl, other_centers(i, :));
   
  % if d_2other < d_2prime
        lb = prob.xl_bl;
        ub = prob.xl_bu;
        
        % fprintf('problems: %s seed %d has problem\n', prob.name, seed_id);
        
        if correctionFunc_decision
            plot2d_withSurrogate(xu, lb, ub, match_xl, seed_xl, prob, seed_xu, delta_mdl, close_mdl, corr_xl);
        else
            % plot2dlower(xu, lb, ub, match_xl, seed_xl, prob, seed_xu);         
            % plot2d_withSurrogate(xu, lb, ub, match_xl, corr_xl, prob, seed_xu, delta_mdl, close_mdl, corr_xl);      
            plot2d_withCokring(xu, lb, ub, match_xl, seed_xl, prob, seed_xu, cokrg_mdl, cokrg_xl, trainxl, trainfl, ...
                co_cheapx,co_cheapy, co_expensivex, co_expensivey);
        end
        
       
 % end
% end

end