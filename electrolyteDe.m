%% Electrolyte Diffusion Coefficient Function: D_e(c_e) [m^2/s]
%   Created July 12, 2011 by Scott Moura

function [D_e,varargout] = electrolyteDe(c_e)

% From DUALFOIL LiPF6 in EC:DMC, Capiglia et al. 1999
D_e = 5.34e-10*exp(-0.65*c_e/1e3);

if(nargout == 2)
    dD_e = -0.65*D_e/1e3;
    varargout{1} = dD_e;
end