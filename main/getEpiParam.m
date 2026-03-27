% Specify an infectious disease and compute its parameters
function [w0, R0, wstat] = getEpiParam(epiNo)

% Features and assumptions
% - returns generation time distribution and R0 for disease
% - generation time w0 is negative binomial
% - remove 1-day shift

% Disease-specific parameters
switch epiNo
    case 1
        % Marburg: van Kerkhove 2015
        epiNam = 'Marburg'; 
        R0 = 3; g0 = 9-1; s = 5.4;
    case 2
        % Measles: Klinkenberg 2011
        epiNam = 'Measles'; R0 = 12; 
        g0 = 12.2-1; s = 3.62;
    case 3
        % COVID-19: Ferguson 2020
        epiNam = 'COVID-19'; R0 = 3.5; 
        g0 = 6.5-1; s = 0.62*6.5;
    case 4
        % Ebola virus: van Kerkhove 2015
        epiNam = 'Ebola virus'; 
        R0 = 2; g0 = 15.3-1; s = 9.3;
    otherwise
        error('Unknown epidemic ID');
end
disp(['Simulating under: ', epiNam, ' parameters']);

% Calculate NB shape parameter
m = (g0^2)/(s^2 - g0); 
% Obtain NB probability parameter
p = m/(g0 + m);  
% Domain of generation time
tw = 0:1000;

% Get NB generation time distribution 
w0 = nbinpdf(tw, m, p); 
% Even though tw(1) = 0, by shift this is w(1 day)

% Statistics of generation time
wm = w0*tw'; wv = w0*(tw.^2)' - wm^2;
% Mean and std deviation
wstat = [wm sqrt(wv)];

% Validate variance matches theory
[~, wv0] = nbinstat(m, p);
if abs(wv - wv0) > 1e-8
    error('Generation time variance mismatch');
end

