%% Electrolyte Activity Coefficient Function: dlnf/dln(c_e) 
%   Created Oct 25, 2016 by Saehong Park

function [dActivity,varargout] = electrolyteAct(c_e,T,p)

% From LiPF6, Valoen et al. 2005
% Fig.6 in the paper

% DataFitting Coefficients

v00 = 0.601;
v01 = 0;
v10 = -0.24;
v11 = 0;
v20 = 0;
v21 = 0;
v30 = 0.982;
v31 = -0.0052;

c_e = c_e/1000; % UnitConversion: 1 mol/L -> 1000 mol/m^3

dActivity = ((v00 + v10.*(c_e).^(0.5) + v30*(1+v31*(T - p.T_ref)) .* (c_e).^(1.5))./(1-p.t_plus))-1;

if(nargout == 2)
    d_dactivity = (0.5 * v10 * (c_e).^(-0.5) + (1.5)*v30*(1+v31*(T - p.T_ref))*(c_e).^(0.5))/(1-p.t_plus);
    varargout{1} = d_dactivity;
end
