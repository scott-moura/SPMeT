%% Li Diffusion in Electrolyte Phase, c_e(x,t) for SPMe
%   Created June 21, 2015 by Scott Moura
%   Adopted from DFN model

function c_ex = electrolyte_scott(data,p)
    
    % run simulation
    x0 = data.c_e0;
    t = data.time;
    [t,x] = ode23s(@(t,x) ode_electrolyte(t,x,data,p),t,x0);
    
    % Parse output
    c_e = x';
    
    c_ex = zeros(p.Nx+1,length(t));
    
    for k = 1:length(t)
    
        % Compute Boundary Conditions
        c_e_bcs = p.ce.C * c_e(:,k);

        % Separate and aggregate
        c_en = c_e(1:(p.Nxn-1),k);
        c_es = c_e((p.Nxn-1)+1:(p.Nxn-1)+(p.Nxs-1),k);
        c_ep = c_e((p.Nxn-1)+p.Nxs : end,k);
        c_ex(:,k) = [c_e_bcs(1); c_en; c_e_bcs(2); c_es; c_e_bcs(3); c_ep; c_e_bcs(4)];
        
    end

%     % run simulation
%     c0 = data.c_e0;
%     c_e_nobd = ode(@ode_electrolyte,data.time,c0,odeset,data,p);
%     
%     % add in boundary condition points
%     h_n = p.L_n / n(1);     % space step size in anode 
%     h_s = p.L_s / n(2);     % space step size in separator
%     h_p = p.L_p / n(3);     % space step size in cathode
%     k_as = (p.epsilon_e_s*h_n)/(p.epsilon_e_n*h_s);
%     k_sc = (p.epsilon_e_p*h_s)/(p.epsilon_e_s*h_p);
%     
%     bcas = inv(1 + k_as)*(c_e_nobd(n(1)-2,:) + k_as*c_e_nobd(n(1)-1));
%     bcsc = inv(1 + k_sc)*(c_e_nobd(n(1)+n(2)-4,:) + k_sc*c_e_nobd(n(1)+n(2)-3,:));
%     
%     c_e = [c_e_nobd(1,:); 
%             c_e_nobd(1:n(1)-2,:); 
%             bcas; 
%             c_e_nobd((n(1)-1):(n(1)+n(2)-4),:);
%             bcsc;
%             c_e_nobd((n(1)+n(2)-3):end,:)
%             c_e_nobd(end,:)];

end

function [c_e_dot] = ode_electrolyte(t,c_e,data,p)

    % Parse and interpolate input current, compute j_n
    cur = lininterp1f(data.time,data.cur,t,[]);
    jn = cur/(p.Faraday*p.a_s_n*p.Area*p.L_n);
    jp = -cur/(p.Faraday*p.a_s_p*p.Area*p.L_p);

    % Compute Boundary Conditions
    c_e_bcs = p.ce.C * c_e;

    % Separate and aggregate
    c_en = c_e(1:(p.Nxn-1));
    c_es = c_e((p.Nxn-1)+1:(p.Nxn-1)+(p.Nxs-1));
    c_ep = c_e((p.Nxn-1)+p.Nxs : end);
    c_ex = [c_e_bcs(1); c_en; c_e_bcs(2); c_es; c_e_bcs(3); c_ep; c_e_bcs(4)];

    % Compute Electrolyte Diffusion Coefficient and Derivative
    [D_en,dD_en] = electrolyteDe(c_en);
    [D_es,dD_es] = electrolyteDe(c_es);
    [D_ep,dD_ep] = electrolyteDe(c_ep);

    % System Matrices have all been precomputed & stored in param struct "p"
    
    % Compute derivative
    c_en_dot = dD_en.*(p.ce.M1n*c_en + p.ce.M2n*c_e_bcs(1:2)).^2 ...
        + D_en.*(p.ce.M3n*c_en + p.ce.M4n*c_e_bcs(1:2)) + diag(p.ce.M5n)*jn;

    c_es_dot = dD_es.*(p.ce.M1s*c_es + p.ce.M2s*c_e_bcs(2:3)).^2 ...
        + D_es.*(p.ce.M3s*c_es + p.ce.M4s*c_e_bcs(2:3));

    c_ep_dot = dD_ep.*(p.ce.M1p*c_ep + p.ce.M2p*c_e_bcs(3:4)).^2 ...
        + D_ep.*(p.ce.M3p*c_ep + p.ce.M4p*c_e_bcs(3:4)) + diag(p.ce.M5p)*jp;

    % Assemble c_e_dot
    c_e_dot = [c_en_dot; c_es_dot; c_ep_dot];
    

%     % n = [n_n,n_s,n_p], is number of discretization points in space
% %     disp(t);
%     if(mod(t,10) == 0)
%         fprintf('.');
%     end
%     cur = data.cur(find(t >= data.time,1,'last'));
%     
%     % compute discretization parameters
%     h_n = p.L_n / n(1);     % space step size in anode 
%     h_s = p.L_s / n(2);     % space step size in separator
%     h_p = p.L_p / n(3);     % space step size in cathode
%     
%     c_e_dot = zeros(sum(n)-6,1);
%     
%     % computing c_e_dot for interior of anode (n)
%     for i=2:n(1)-3
%         c_e_dot(i) = (p.D_e/h_n^2)*(c_e(i-1) - 2*c_e(i) + c_e(i+1)) + ...
%             ((1-p.t_plus)/(p.epsilon_e_n*p.Faraday*p.Area*p.L_n))*cur;
%     end
%     % computing c_e_dot for interior of separator (s)
%     for i=(2:n(2)-3)+(n(1)-2)
%         c_e_dot(i) = (p.D_e/h_s^2)*(c_e(i-1) - 2*c_e(i) + c_e(i+1));
%     end
%     % computing c_e_dot for interior of cathode (p)
%     for i=(2:n(3)-3)+(n(1)-2)+(n(2)-2)
%         c_e_dot(i) = (p.D_e/h_p^2)*(c_e(i-1) - 2*c_e(i) + c_e(i+1)) - ...
%             ((1-p.t_plus)/(p.epsilon_e_p*p.Faraday*p.Area*p.L_p))*cur;
%     end
%     
%     % computing boundary cases:
%     % anode-0 (0-): zero flux
%     c_e_dot(1) = (p.D_e/h_n^2)*(c_e(2) - c_e(1)) + ... 
%         ((1-p.t_plus)/(p.epsilon_e_n*p.Faraday*p.Area*p.L_n))*cur;
%     % cathode-0 (0+): zero flux
%     c_e_dot(end) = (p.D_e/h_p^2)*(-c_e(end) + c_e(end-1)) - ...
%         ((1-p.t_plus)/(p.epsilon_e_p*p.Faraday*p.Area*p.L_p))*cur;
%     % anode/separator boundaries (L-, 0sep): continuity
%     k_as = (p.epsilon_e_s*h_n)/(p.epsilon_e_n*h_s);
%     c_e_dot(n(1)-2) = (p.D_e/h_n^2)*(c_e(n(1)-3) - 2*c_e(n(1)-2) + ...
%         inv(1+k_as)*(c_e(n(1)-2) + k_as*c_e(n(1)-1))) + ...
%         ((1-p.t_plus)/(p.epsilon_e_n*p.Faraday*p.Area*p.L_n))*cur;
%     c_e_dot(n(1)-1) = (p.D_e/h_s^2)*(c_e(n(1)) - 2*c_e(n(1)-1) + ...
%         inv(1+k_as)*(c_e(n(1)-2) + k_as*c_e(n(1)-1)));
%     % cathode/separator boundaries (L+, Lsep): continuity
%     k_sc = (p.epsilon_e_p*h_s)/(p.epsilon_e_s*h_p);
%     c_e_dot(n(1)+n(2)-4) = (p.D_e/h_s^2)*(c_e(n(1)+n(2)-5) - 2*c_e(n(1)+n(2)-4) + ...
%         inv(1+k_sc)*(c_e(n(1)+n(2)-4) + k_sc*c_e(n(1)+n(2)-3)));
%     c_e_dot(n(1)+n(2)-3) = (p.D_e/h_p^2)*(c_e(n(1)+n(2)-2) - 2*c_e(n(1)+n(2)-3) + ...
%         inv(1+k_sc)*(c_e(n(1)+n(2)-4) + k_sc*c_e(n(1)+n(2)-3))) - ...
%         ((1-p.t_plus)/(p.epsilon_e_p*p.Faraday*p.Area*p.L_p))*cur;
end