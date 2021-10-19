%-----------------------
% debug parameter
visual          = false;
coresteps       = false; % true means no upper local search or re-evaluation
%------------------------
% outer loop parameter
inisize_u         = 20;   % upper level initialization sample size
inisize_l         = 30;   % lower level initialization sample size
numiter_l         = 70;   % lower level infill iteration number
numiter_u         = 20;   % upper level infill iteration number
num_pop           = 100;  % EA search on surrogate( KE or EI), population size
num_gen           = 100;  % EA search on surrogate( KE or EI), generation size
maxFE_l           = 200;
%---------------------------
% bilevel local search parameter
blsearch                     =  false; % post infill process: whether conduct bilevel local search on upper level
nest_parameter.maxUpperFE    =  100; % post infill process: upper level local search (SQP) allowed number of FEs
nest_parameter.maxLowerFE    =  100; % post infill process: nest lower level local search(SQP) allowed number of FEs
nest_parameter.global_search =  true;% post infill process: bilevel local search on the lower level EA(global_search) + SQP
nest_parameter.num_gen       =  5;  % post infill process: if lower level EA + SQP, EA's generation
nest_parameter.num_pop       =  10;   % post infill process: if lower level EA + SQP, EA's population size

%----------------------------------
% re-evaluation parameter
re_evalparameter.maxLowerFE    =  100;   % post infill process: re-evaluation on lowerlevel local search (SQP), max allowed FEs 
re_evalparameter.global_search =  true;  % post infill process: re-evaluation EA + SQP
re_evalparameter.num_gen       =  5;    % post infill process: re-evaluation EA generation
re_evalparameter.num_pop       =  10;     % post infill process: re-evaluation EA population

%------------------------------------
llmatch_p               = struct();      % llmatch function parameter
llmatch_p.num_pop       = num_pop;
llmatch_p.num_gen       = num_gen;
llmatch_p.egoinitsize   = inisize_l;
llmatch_p.egoitersize   = numiter_l;
llmatch_p.maxFE_l      = maxFE_l;
