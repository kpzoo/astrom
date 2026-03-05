% Initialisation settings for plotting and saving
function [gr1, gr2, cmap, cols] = initialisePlot(nCol, fnt)

% Features and assumptions
% - adjusts figure font sizes


% Ensure linspecer included
if exist('linspecer', 'file') ~= 2
    error('The function linspecer is not found. Add it to your path.');
end

% Set global graphics defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter',...
    'latex', 'defaultTextInterpreter', 'latex', 'defaultAxesFontSize', fnt,...
    'defaultAxesYGrid', 'off', 'defaultAxesXGrid', 'off');

% Colour map
cmap = linspecer(nCol, 'sequential'); 
% Grey colours
gr1 = 0.8*ones(1, 3); gr2 = 0.5*ones(1, 3);
% Default color cell array
cols = {'b', gr1, gr2, 'g', 'm', 'r'};