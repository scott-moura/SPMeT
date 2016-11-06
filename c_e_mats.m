%% Matrices for Li Diffusion in Electrolyte Phase, c_e(x,t)
%   Created June 19, 2015 by Scott Moura
%   A complete reboot of the original c_e_mats from July 2011

function [M1n,M2n,M3n,M4n,M5n, M1s,M2s,M3s,M4s, M1p,M2p,M3p,M4p,M5p, C_ce] = c_e_mats(p)

%% Lumped Coefficients
Del_xn = p.L_n * p.delta_x_n;
Del_xs = p.L_s * p.delta_x_s;
Del_xp = p.L_p * p.delta_x_p;

%% Matrices in nonlinear dynamics
M1n = sparse((diag(ones(p.Nxn-2,1),+1) - diag(ones(p.Nxn-2,1),-1))/(2*Del_xn));
M1s = sparse((diag(ones(p.Nxs-2,1),+1) - diag(ones(p.Nxs-2,1),-1))/(2*Del_xs));
M1p = sparse((diag(ones(p.Nxp-2,1),+1) - diag(ones(p.Nxp-2,1),-1))/(2*Del_xp));


M2n = zeros(p.Nxn-1,2);
M2n(1,1) = -1/(2*Del_xn);
M2n(end,end) = 1/(2*Del_xn);
M2n = sparse(M2n);

M2s = zeros(p.Nxs-1,2);
M2s(1,1) = -1/(2*Del_xs);
M2s(end,end) = 1/(2*Del_xs);
M2s = sparse(M2s);

M2p = zeros(p.Nxp-1,2);
M2p(1,1) = -1/(2*Del_xp);
M2p(end,end) = 1/(2*Del_xp);
M2p = sparse(M2p);


M3n = sparse((-2*diag(ones(p.Nxn-1,1),0) + diag(ones(p.Nxn-2,1),+1) + diag(ones(p.Nxn-2,1),-1))/(Del_xn^2));
M3s = sparse((-2*diag(ones(p.Nxs-1,1),0) + diag(ones(p.Nxs-2,1),+1) + diag(ones(p.Nxs-2,1),-1))/(Del_xs^2));
M3p = sparse((-2*diag(ones(p.Nxp-1,1),0) + diag(ones(p.Nxp-2,1),+1) + diag(ones(p.Nxp-2,1),-1))/(Del_xp^2));


M4n = zeros(p.Nxn-1,2);
M4n(1,1) = 1/(Del_xn^2);
M4n(end,end) = 1/(Del_xn^2);
M4n = sparse(M4n);

M4s = zeros(p.Nxs-1,2);
M4s(1,1) = 1/(Del_xs^2);
M4s(end,end) = 1/(Del_xs^2);
M4s = sparse(M4s);

M4p = zeros(p.Nxp-1,2);
M4p(1,1) = 1/(Del_xp^2);
M4p(end,end) = 1/(Del_xp^2);
M4p = sparse(M4p);

M5n = (1-p.t_plus)*p.a_s_n/p.epsilon_e_n * speye(p.Nxn-1);
M5p = (1-p.t_plus)*p.a_s_p/p.epsilon_e_p * speye(p.Nxp-1);

%% Boundary Conditions
N1 = zeros(4,p.Nx-3);
N2 = zeros(4);

% BC1
N1(1,1) = +4;
N1(1,2) = -1;
N2(1,1) = -3;

% BC2
N1(2,p.Nxn-2) = (p.epsilon_e_n^p.brug)/(2*Del_xn);
N1(2,p.Nxn-1) = (-4*p.epsilon_e_n^p.brug)/(2*Del_xn);
N2(2,2) = (3*p.epsilon_e_n^p.brug)/(2*Del_xn) + (3*p.epsilon_e_s^p.brug)/(2*Del_xs);
N1(2,p.Nxn) = (-4*p.epsilon_e_s^p.brug)/(2*Del_xs);
N1(2,p.Nxn+1) = (p.epsilon_e_s^p.brug)/(2*Del_xs);

% BC3
N1(3,p.Nxn+p.Nxs-3) = (p.epsilon_e_s^p.brug)/(2*Del_xs);
N1(3,p.Nxn+p.Nxs-2) = (-4*p.epsilon_e_s^p.brug)/(2*Del_xs);
N2(3,3) = (3*p.epsilon_e_s^p.brug)/(2*Del_xs) + (3*p.epsilon_e_p^p.brug)/(2*Del_xp);
N1(3,p.Nxn+p.Nxs-1) = (-4*p.epsilon_e_p^p.brug)/(2*Del_xp);
N1(3,p.Nxn+p.Nxs) = (p.epsilon_e_p^p.brug)/(2*Del_xp);


% BC4
N1(4,end-1) = +1;
N1(4,end) = -4;
N2(4,4) = +3;

%%% SPARSE OUTPUT 
C_ce = sparse(-N2\N1);


