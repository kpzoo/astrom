%% Perturb generation time distribution into a new one
function [w1, R1, wtest] = perturbW(w0, R0, ndim, r0, wdiff)

% Features and assumptions
% - constrain mean of w, growth rate r and reproduction number R
% - change w0 to w1 with cost wdiff under constraints

% Structured population matrices
T = diag(ones(1, ndim-1), -1); 
B = [1 zeros(1, ndim-1)]'; In = eye(ndim);

% Bounds on the new distribution probabilities
lb = zeros(1, ndim); ub = ones(1, ndim);
% Ensure tail probabilities not too large
maxDay = 31; ub(maxDay+1:maxDay+7) = 3*w0(maxDay);
ub(maxDay+7:ndim) = 0.001;

% Also works with maxDay = 40

% Equality constraints on mean and sum of w
Aeq = [ones(1, ndim); 1:ndim];
beq = [1; sum(w0.*(1:ndim))];

% Initial uniform distribution and set options
init = ones(1, ndim)/ndim;
fmopts = optimoptions('fmincon', 'Algorithm','sqp', 'StepTolerance', 1e-12,...
    'ConstraintTolerance', 1e-10, 'MaxFunctionEvaluations', 2e5);

% Run optimisation to get new generation time distribution
w1 = fmincon(@(w) minW(w, w0, R0, R0, B, T, In, r0, wdiff), init, [], [],...
    Aeq, beq, lb, ub, [], fmopts);

% Properties of optimal solution
A = R0*w1; F = B*A; r1 = max(eig(F + T)); R1 = R0;
wtest.r = [r0 r1]; wtest.mean = (1:ndim)*[w0' w1']-1;


%% Objective function to minimise
function cost = minW(w, w0, R, R0, B, T, In, r0, wdiff)

% Construct fertility matrix
A = R*w; F = B*A;
% Growth rate 
r = max(eig(F + T));
% NGM or reproduction number
NGM = A*((In - T)\B);

% Cost to keep r close to r0
cost = abs(r - r0)*100;
% Cost to keep R close to R0
cost = cost + abs(NGM - R0);
% But cannot be too close to w0 (KL divergence)
cost = cost + wdiff*sum(w .* log((w+1e-12)./(w0+1e-12)));