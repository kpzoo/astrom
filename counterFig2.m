%% Two epidemics with OL divergence but CL convergence
clearvars; close all; clc; 

% Features and assumptions
% - reproduces Fig 2 of paper
% - control attenuates structural uncertainties
% - parameters from Ebola virus and COVID-19

% Dependencies and general code
mainDir = fullfile(pwd, 'main'); addpath(genpath(mainDir));
% MPC specific codebase (not toolbox)
mpcDir = fullfile(pwd, 'mpc'); addpath(genpath(mpcDir));
% Default plotting options
fnt = 22; [gr1, gr2, cmap, cols] = initialisePlot(20, fnt);

%% Simulate different diseases OL and CL

% Possible diseases to simulate
epiNam = {"Marburg", "Measles", "COVID", "Ebola"}; 
epiID = 3:4; lenID = length(epiID);

% Plotting and if want contiguous control type
indivPlt = 0; contigtype = 1;

% Main loop across diseases
dis = cell(1, lenID);
for ii = 1:lenID
    % Trajectories and statistics for disease
    dis{ii} = getCLSynch(epiID(ii), fnt, indivPlt, contigtype);
    dis{ii}.name = epiNam{epiID(ii)};
end

%% Plot results for both diseases

figure('Color','w','Position',[160 50 1000 1000]);
inset0 = [0.32 0.68 0.1 0.08];

% Incidence trajectories for OL and then CL
for ii = 1:2
    sim = dis{ii};
    subplot(4, 2, [ii ii+2]);
    stairs(sim.t, sim.I0, 'b', 'LineWidth', 2);
    hold on; grid on;
    stairs(sim.t, sim.I1, 'r', 'LineWidth', 2);
    h = gca; yl = h.YLim;
    plot([sim.tst sim.tst], yl, 'k--', 'LineWidth', 2);
    hold off; box off; title(dis{ii}.name, FontSize=fnt);
    xlabel('$t$ (days)', 'FontSize',fnt);
    ylabel('$i(t)$', 'FontSize',fnt);

    % Inset axes of cumulative infectons
    switch ii
        case 1
            axes('Position', inset0);
        case 2
            inset1 = inset0; inset1(1) = inset1(1) + 0.45;
            axes('Position', inset1);
    end
    plot(sim.t, cumsum(sim.I0), 'b-', 'LineWidth', 2);
    hold on; grid on;
    plot(sim.t, cumsum(sim.I1), 'r-', 'LineWidth', 2);
    h = gca; yl = h.YLim;
    plot([sim.tst sim.tst], yl, 'k--', 'LineWidth', 2);
    hold off; box off;
    title('$\sum_t i(t)$', FontSize = fnt);
end

% Generation time distributions and control
inset2 = [0.32 0.38 0.1 0.08];
for ii = 1:2
    sim = dis{ii};
    subplot(4, 2, ii+4);
    stairs(sim.tdim-1, sim.w0, 'b-', 'LineWidth', 2);
    hold on; box off;
    stairs(sim.tdim-1, sim.w1, 'r-', 'LineWidth', 2);
    grid on; hold off; 
    xlim([0 sim.tdim(end)]);
    xlabel('$s$ (days)', 'FontSize',fnt);
    ylabel('$w(s)$', 'FontSize',fnt);

    if max(abs(sim.g0 - sim.g1)) < 10^-6
        hold on; h = gca; yl = h.YLim;
        plot([sim.g0 sim.g0], yl, 'k--', 'LineWidth', 2);
        hold off; 
    else
        % The mean generation times differ
        hold on; h = gca; yl = h.YLim;
        plot([sim.g0 sim.g0], yl, 'b--', 'LineWidth', 2);
        plot([sim.g1 sim.g1], yl, 'r--', 'LineWidth', 2);
        hold off;
    end

    % Inset axes of control proportions removed
    switch ii
        case 1
            axes('Position', inset2);
        case 2
            inset3 = inset2; inset3(1) = inset3(1) + 0.45;
            axes('Position', inset3);
    end
    stairs(sim.tdim-1, sim.e, 'k-', 'LineWidth', 1);
    grid on; box off; xlim([0 sim.tdim(end)]);
    title('control', FontSize = fnt);
end

% Statistics of growth and reproduction
for ii = 1:2
    sim = dis{ii};
    subplot(4, 2, ii+6);
    hold on; grid off; msize = 6;

    % Shade subcritical region (R < 1, r < 1)
    xr = xlim; yr = ylim;

    patch([xr(1) 1 1 xr(1)], [yr(1) yr(1) 1 1],[0.9 0.9 0.9], ...
        'EdgeColor','none','FaceAlpha',0.8);

    plot(sim.R0, sim.r0, 'bs', 'LineWidth', 2, 'MarkerSize', msize, 'MarkerFaceColor', 'b');
    plot(sim.R1, sim.r1, 'rs', 'LineWidth', 2, 'MarkerSize', msize, 'MarkerFaceColor', 'r');
    plot(sim.Rk0, sim.rk0, 'bo', 'LineWidth', 2, 'MarkerSize', msize, 'MarkerFaceColor', 'b');
    plot(sim.Rk1, sim.rk1, 'ro', 'LineWidth', 2, 'MarkerSize', msize, 'MarkerFaceColor', 'r');

    plot([sim.R0 sim.Rk0],[sim.r0 sim.rk0], 'b-', 'LineWidth', 1);
    plot([sim.R1 sim.Rk1],[sim.r1 sim.rk1], 'r-', 'LineWidth', 1);
    xline(1, 'k--', 'LineWidth', 1);
    yline(1, 'k--', 'LineWidth', 1);

    xlabel('$R$', 'FontSize', fnt);
    ylabel('$r$', 'FontSize', fnt);
    box off; grid on; hold off;

    xlim([min(sim.Rk0, sim.Rk1), max(sim.R0, sim.R1)]);
    ylim([min(sim.rk0, sim.rk1), max(sim.r0, sim.r1)]);
end