%% Li Diffusion in Electrolyte Phase, c_e(x,t) for SPMe
%   Created June 21, 2015 by Scott Moura
%   Adopted from DFN model

function c_ex = electrolyte_scott(data,p)
    
    Nn = p.Nxn - 1;
    Np = p.Nxp - 1;
    Nx = p.Nx - 3;

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
        c_en = c_e(1:Nn,k);
        c_es = c_e(Nn+1:Nn+p.Nxs-1,k);
        c_ep = c_e(Nn+p.Nxs : end,k);
        c_ex(:,k) = [c_e_bcs(1); c_en; c_e_bcs(2); c_es; c_e_bcs(3); c_ep; c_e_bcs(4)];
        
    end

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

    % ADD BRUGGEMAN RELATION % Apr.22 2016 by Saehong Park
    D_en_eff = D_en .* p.epsilon_e_n.^(p.brug-1);
    dD_en_eff = dD_en .* p.epsilon_e_n.^(p.brug-1);

    % DO same s,p 
    D_es_eff = D_es .* p.epsilon_e_s.^(p.brug-1);
    dD_es_eff = dD_es .* p.epsilon_e_s.^(p.brug-1);

    D_ep_eff = D_ep .* p.epsilon_e_p.^(p.brug-1);
    dD_ep_eff = dD_ep .* p.epsilon_e_p.^(p.brug-1);
    
    % System Matrices have all been precomputed & stored in param struct "p"
    
    % Compute derivative
    c_en_dot = dD_en_eff.*(p.ce.M1n*c_en + p.ce.M2n*c_e_bcs(1:2)).^2 ...
        + D_en_eff.*(p.ce.M3n*c_en + p.ce.M4n*c_e_bcs(1:2)) + diag(p.ce.M5n)*jn;

    c_es_dot = dD_es_eff.*(p.ce.M1s*c_es + p.ce.M2s*c_e_bcs(2:3)).^2 ...
        + D_es_eff.*(p.ce.M3s*c_es + p.ce.M4s*c_e_bcs(2:3));

    c_ep_dot = dD_ep_eff.*(p.ce.M1p*c_ep + p.ce.M2p*c_e_bcs(3:4)).^2 ...
        + D_ep_eff.*(p.ce.M3p*c_ep + p.ce.M4p*c_e_bcs(3:4)) + diag(p.ce.M5p)*jp;

    % Assemble c_e_dot
    c_e_dot = [c_en_dot; c_es_dot; c_ep_dot];
    

end