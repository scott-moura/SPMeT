%% Compute Initial Solid Concentrations from Voltage and Params
%   Created Mar 12, 2014

function [csn0,csp0] = init_cs(p,V)

%% Use Bisection Algorithm

% Algorithm params
maxiters = 50;
x = zeros(maxiters,1);
f = nan*ones(maxiters,1);
tol = 1e-5;

% Initial Guesses
x_low = 0.2 * p.c_s_p_max;
x_high = 1.0 * p.c_s_p_max;
x(1) = 0.6 * p.c_s_p_max;

% Iterate Bisection Algorithm
for idx = 1:maxiters

    theta_p = x(idx)/p.c_s_p_max;
    theta_n = (p.n_Li_s-p.epsilon_s_p*p.L_p*p.Area*x(idx))/(p.c_s_n_max*p.epsilon_s_n*p.L_n*p.Area);

    OCPn = refPotentialAnode(p,theta_n);
    OCPp = refPotentialCathode(p,theta_p);

    f(idx) = OCPp - OCPn - V;
        
    if(abs(f(idx)) <= tol)
        break;
    elseif(f(idx) <= 0)
        x_high = x(idx);
    else
        x_low = x(idx);
    end

    % Bisection
    x(idx+1) = (x_high + x_low)/2;
    x(idx+1)/p.c_s_p_max;

end

% Output conveged csp0
csp0 = x(idx);

%% Use Newton-Raphson to Solve OCP functions
% 
% % Initial Guess for csp0
% x0 = 0.6 * p.c_s_p_max;
% 
% % Iterate!
% max_iters = 50;
% tol = 1e-5;
% 
% % Preallocate
% % x = nan*ones(max_iters,1);
% % x(1) = x0;
% % f = zeros(max_iters,1);

% for idx = 1:max_iters
% 
%     theta_p = x(idx)/p.c_s_p_max;
%     theta_n = (p.n_Li_s-p.epsilon_s_p*p.L_p*x(idx))/(p.c_s_n_max*p.epsilon_s_n*p.L_n);
% 
%     [OCPn,dOCPn] = refPotentialAnode(p,theta_n);
%     [OCPp,dOCPp] = refPotentialCathode(p,theta_p);
% 
%     f(idx) = OCPp - OCPn - V;
%     df = dOCPp/p.c_s_p_max + (p.epsilon_s_p*p.L_p)/(p.c_s_n_max*p.epsilon_s_n*p.L_n) * dOCPn;
% 
%     Del_x = -f(idx)/df;
% 
%     x(idx+1) = x(idx) + Del_x;
%     
% %     if(f <= tol)
% %         break;
% %     end
%     
% end
% 
% x/p.c_s_p_max
% f
% 
% % Output conveged csp0
% csp0 = x(idx);

%% Compute csn0
csn0 = (p.n_Li_s - p.epsilon_s_p * p.L_p * p.Area * csp0) / (p.epsilon_s_n * p.L_n * p.Area);