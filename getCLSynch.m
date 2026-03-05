%% Generate two epidemics with OL convergence but CL divergence
function sim = getCLSynch(epiID, fnt, indivPlt, contigtype)

% Features and assumptions
% - simulates two structured renewal models with control
% - control by removing infections in contiguous ages
% - models have similar generation time distributions

%% Shared epidemic parameter settings

% Possible diseases to simulate
epiNam = {"Marburg", "Measles", "COVID", "Ebola"};
% Choose a disease and get parameters
[w0, R0, ~] = getEpiParam(epiID);

% Settings for showing control
if strcmp(epiNam{epiID}, 'Ebola')
    
    % Time settings
    ndim = 80; tst = 120; T = 300;
    % Optimisation settings
    thresh = 0.016; wdiff = 0.0021;

elseif strcmp(epiNam{epiID}, 'COVID')
    
    % Time settings
    ndim = 60; tst = 40; T = 150; 
    % Optimisation settings
    thresh = 0.03; wdiff = 0.003;

else
    error('Disease not yet supported');
end
    
% Time series and initialise infections
tdim = 1:ndim; t = 1:T; Iinit = 2;

% Baseline epidemic system
[X0, Q0, w0, ~, F0, A0] = initStateMx(tdim, ndim, T, Iinit, R0, w0);
% Population model and poles
pop0 = getPopModel(A0, ndim); r0 = pop0.r(1);

%% Construct second system with different OL

% Alter w0 into w1 (slightly uncertain, reshapes mass)
switch contigtype
    case 0
        % Approach using a Gaussian bump second mode
        if strcmp(epiNam{epiID}, 'Ebola')

            % Contiguous late-age block for second mode
            id1 = find(w0(2:end) < 0.1, 1, 'first') + 1;
            id1 = min(id1, ndim-15); id2 = id1 + 14;

            % Region will add infectiousness mass eps
            idadd = id1:id2; eps = 0.2;
            % Shape a Gaussian around the idadd region
            mu = round(mean(idadd)); sig = (id2 - id1)/3;

            % Gaussian bump in mass
            bump = exp(-(tdim - mu).^2/(2*sig^2));
            bump = bump.*ismember(tdim, idadd);
            bump = bump/sum(bump); w1 = w0 + eps*bump;

        elseif strcmp(epiNam{epiID}, 'COVID')

            % Contiguous late-age block for second mode
            id1 = find(w0(2:end) < 0.1, 1, 'first') + 1;
            id1 = min(id1, ndim-15); id2 = id1 + 14;

            % Region will add infectiousness mass eps
            idadd = id1:id2; eps = 0.25;
            % Shape a Gaussian around the idadd region
            mu = round(mean(idadd)); sig = (id2 - id1)/3;

            % Gaussian bump in mass
            bump = exp(-(tdim - mu).^2/(2*sig^2));
            bump = bump.*ismember(tdim, idadd);
            bump = bump/sum(bump); w1 = w0 + eps*bump;

        end

    case 1
        % Approach using different bump types
        if strcmp(epiNam{epiID}, 'Ebola')

            % Contiguous late-age block for second mode
            id1 = find(w0 >= 0.02, 1, 'first');
            id1 = min(id1, ndim-15); id2 = id1 + 14;

            % Region will add infectiousness mass eps
            idadd = id1:id2; eps = 0.1;
            % Shape a Gaussian around the idadd region
            mu = round(mean(idadd)); sig = (id2 - id1)/3;

            % Gaussian bump in mass
            bump = exp(-(tdim - mu).^2/(2*sig^2));
            bump = bump.*ismember(tdim, idadd);
            bump = bump/sum(bump); w1 = w0 + eps*bump;

        elseif strcmp(epiNam{epiID}, 'COVID')

            % Contiguous late-age block for second mode
            id1 = find(w0 >= 0.02, 1, 'first');
            id1 = min(id1, ndim-15); id2 = id1 + 14;

            % Region will add infectiousness mass eps
            idadd = id1:id2; eps = 0.2;
            % Shape a Gaussian around the idadd region
            mu = round(mean(idadd)); sig = (id2 - id1)/3;

            % Gaussian bump in mass
            bump = exp(-(tdim - mu).^2/(2*sig^2));
            bump = bump.*ismember(tdim, idadd);
            bump = bump/sum(bump); w1 = w0 + eps*bump;

        end

end

if any(w1 < 0)
    % Test that optimised w1 is feasible
    if abs(min(w1)) < 10^(-20)
        w1(w1 < 0) = 0; w1 = w1/sum(w1);
    else
        error('Check optimised w1 for negative values');
    end
end
    
%% Obtain second OL system and its statistics

% Re-normalise generation time and get R1
A1 = R0*w1; R1 = sum(A1); w1 = w1/sum(w1);
    
% Second population model and poles
[X1, Q1, w1, Tsh, F1, A1] = initStateMx(tdim, ndim, T, Iinit, R1, w1);
pop1 = getPopModel(A1, ndim); r1 = pop1.r(1);

% Simulate uncontrolled (open loop) epidemic until tst
for ii = 2:tst
    % State is all infectious indviduals of different ages
    [X0(:, ii), Q0(:, ii)] = renewalStep(Q0(:, ii-1), X0(:, ii-1),...
        F0, zeros(ndim, ndim), Tsh, 0);
    [X1(:, ii), Q1(:, ii)] = renewalStep(Q1(:, ii-1), X1(:, ii-1),...
        F1, zeros(ndim, ndim), Tsh, 0);
end


%% Compute characteristics of control projections

% Proportion removed (1-1/R0 is critical)
e = (1-1/(R0^2))*ones(1, ndim); 

% Strength based on disease
    switch epiID
        case 3
            % COVID removal proportions
            str = 0.99;
        case 4
            % Ebola removal proportions
            str = 0.99;
    end

% In region added mass increase control
e(idadd) = str; e = e(1:ndim); E = diag(e);

%% Apply the same controller to both epidemics

% Resulting row A matrices that are subtracted
Ae0 = A0.*(1-e); Ae1 = A1.*(1-e);

% Closed loop control applied from tst+1 to end T
for ii = tst+1:T
    % Apply control action and obtain new state and removals
    [X0(:, ii), Q0(:, ii)] = renewalStep(Q0(:, ii-1), X0(:, ii-1),...
        F0, E, Tsh, 1);
    [X1(:, ii), Q1(:, ii)] = renewalStep(Q1(:, ii-1), X1(:, ii-1),...
        F1, E, Tsh, 1);
end

% Incidence curves for both epidemics
I0 = X0(1,:); I1 = X1(1,:); 
% Mean generation times for both epidemics
g0 = w0*tdim'-1; g1 = w1*tdim'-1;

% Asymptotic statistics of controlled systems
popk0 = getPopModel(Ae0, ndim); rk0 = popk0.r; Rk0 = popk0.R;
popk1 = getPopModel(Ae1, ndim); rk1 = popk1.r; Rk1 = popk1.R;

%% Store key epidemic data for plotting

% Incidence trajectories
sim.I0 = I0; sim.I1 = I1;
% Timing variables
sim.t = t; sim.tst = tst; sim.tdim = tdim;
% Generation time distributions
sim.w0 = w0; sim.w1 = w1;
% Mean generation times
sim.g0 = g0; sim.g1 = g1;

% OL growth rates
sim.r0 = r0; sim.r1 = r1;
% CL growth rates
sim.rk0 = rk0; sim.rk1 = rk1;
% OL reproduction numbers
sim.R0 = R0; sim.R1 = R1;
% CL reproduction numbers
sim.Rk0 = Rk0; sim.Rk1 = Rk1;
% Control settings and action
sim.e = e; sim.A0 = A0; sim.A1 = A1;
sim.Ae0 = Ae0; sim.Ae1 = Ae1;
sim.thresh = thresh; sim.wdiff = wdiff;


%% Visualise generation times and trajectories

if indivPlt

    figure('Color','w','Position',[200 200 800 800]);
    tiledlayout(4,2, 'TileSpacing','compact', 'Padding','compact');

    % OL and CL incidence trajectories
    nexttile([2 2])
    stairs(t, I0, 'b', 'LineWidth', 2);
    hold on; grid off;
    stairs(t, I1, 'r', 'LineWidth', 2);
    h = gca; yl = h.YLim;
    plot([tst tst], yl, 'k--', 'LineWidth', 2);
    text(tst, yl(2) - 0.01*(yl(2)-yl(1)), '$\tau$',...
        'Interpreter','latex', 'HorizontalAlignment','center',...
        'VerticalAlignment','bottom', 'FontSize', fnt);
    box off; hold off;
    xlabel('$t$ (days)', 'FontSize',fnt);
    ylabel('$i(t)$', 'FontSize',fnt);
    lgd = legend(epiNam{epiID}, 'Approx model', 'Location', 'best', fontsize = fnt);
    set(lgd, 'Box', 'off');

    % Generation times
    nexttile([2 1])
    stairs(tdim, w0, 'b-', 'LineWidth', 2);
    hold on; box off;
    stairs(tdim, w1, 'r-', 'LineWidth', 2);
    h = gca; yl = h.YLim;
    xg = (wtest.mean(1) - 1);
    plot([xg xg], yl, 'k--', 'LineWidth', 2);
    text(xg, yl(2) - 0.01*(yl(2)-yl(1)), '$g$',...
        'Interpreter','latex', 'HorizontalAlignment','center',...
        'VerticalAlignment','bottom', 'FontSize', fnt);
    grid off; hold off;
    xlabel('$s$ (days)', 'FontSize',fnt);
    ylabel('$w(s)$', 'FontSize',fnt);

    % r comparison
    nexttile
    plot([0 1], [r0 rk0], 'x-', 'LineWidth', 2, 'MarkerSize', 10);
    hold on;
    plot([0 1], [r1 rk1], 'x-', 'LineWidth', 2, 'MarkerSize', 10);
    plot([0 1], [1 1], 'k--');
    box off; grid off; hold off;
    set(gca, 'XTick', [0 1], 'XTickLabel', {'Pre','Post'});
    xlabel('$\Delta r$', 'FontSize', fnt);

    % R comparison
    nexttile
    plot([0 1], [R0 Rk0], 'x-', 'LineWidth', 2, 'MarkerSize', 10);
    hold on;
    plot([0 1], [R1 Rk1],'x-', 'LineWidth', 2, 'MarkerSize', 10);
    plot([0 1], [1 1], 'k--');
    box off; grid off; hold off;
    set(gca, 'XTick', [0 1], 'XTickLabel', {'Pre','Post'});
    xlabel('$\Delta R$', 'FontSize', fnt);

end
