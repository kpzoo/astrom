%% Compute uncontrolled state space matrices
function [X, Q, wpop, Tsh, F, Arow] = initStateMx(tdim, ndim, T, totI, R0, wpop)

% Features and assumptions
% - structured population renewal model F, Tsh
% - ndim complexity via generation time distribution wpop

% Truncate generation time, first element is time 0
wsum = sum(wpop(tdim));
if wsum < 0.99
    % Ensure covering enough of generation time mass
    error('Add more initial infection time series');
else
    wpop = wpop(1:ndim); wpop = wpop/wsum;
end

% Top row of state matrix A
Arow = (R0*wpop(:))';
% Shift transition matrix (sparse)
Tsh = spdiags(ones(ndim-1,1), -1, ndim, ndim);
% Fertility matrix so A = F + Tsh
F = sparse(ones(1,ndim), tdim, Arow, ndim, ndim);

% Population and controller states (ages)
X = zeros(ndim, T); Q = X;
% Initialisation of infectious individuals
X(:, 1) = [totI zeros(1, ndim-1)]';
