%% Fast structured population simulation steps
function [Xnext, Qnext] = renewalStep(Qlast, Xlast, F, E, Tsh, hasCtrl)

% Assumptions and notes
% - input uncontrolled population state space vectors 
% - control using diagonal E removing fractions per age
% - output incidence, controlled Q and infectious X population

% Update age-structure only if control on
if hasCtrl
    % State variable for removal actions
    Qnext = Tsh*Qlast; Qnext(1) = Xlast(1);   
else
    % Set to zero without reallocating zeros(ndim,1)
    Qnext = Qlast;
end

% Standard growth of infectious states
Xnext = (F + Tsh)*Xlast;

% Remove infections by age if applying control via E
if hasCtrl
    Xnext = Xnext - F*E*Qnext;
end
