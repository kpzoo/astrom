% Get NGM and other properties from population model
function popMod = getPopModel(A, ndim)

% Assumptions and notes
% - input population incidence matrix A and size ndim

% Shift operator for ages of infections
T = diag(ones(1, ndim-1), -1); 
% Fertility matrix F based on B and A
B = [1 zeros(1, ndim-1)]'; F = B*A;

% Next generation matrix and growth rate
In = eye(ndim); r = max(eig(F + T));
NGM = A*((In - T)\B);

% Store matrices, TFs and metrics
popMod.R = NGM; popMod.r = r;
popMod.A = A; popMod.B = B; 
popMod.F = F; popMod.T = T;
