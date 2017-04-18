%% SPMeT Light and Fast Data
%   Published April 18, 2017 by Professor Scott Moura
%   Energy, Controls, and Applications Lab (eCAL)
%   University of California, Berkeley
%   http://ecal.berkeley.edu/

%   Code based on publications
%   Battery State Estimation for a Single Particle Model with Electrolyte Dynamics 
%   S. J. Moura, F. Bribiesca Argomedo, R. Klein, A. Mirtabatabaei, M. Krstic 
%   IEEE Transactions on Control System Technology, to appear 
%   DOI: 10.1109/TCST.2016.2571663

%   Optimal Charging of Batteries via a Single Particle Model with Electrolyte and Thermal Dynamics 
%   H. Perez, X. Hu, S. J. Moura 
%   2016 American Control Conference
%   DOI: 10.1109/ACC.2016.7525538
fs = 18;


%% Data

% Test name ; Comp time per 1 sec sim time [sec] ; workspace memory per 1 sec sim time [kB]
data  = 
{'1C Dis 2min', 0.0183, 2.26; ...
 '2C Dis 2min', 0.0142, 2.26; ...
 '5C Dis 2min', 0.0142, 2.26; ...
 '1C Chg 2min', 0.0092, 2.26; ...
 '2C Chg 2min', 0.0108, 2.26; ...
 '5C Chg 2min', 0.0142, 2.26; ...
 'UDDSx2 63min', 0.0131, 4.21};
 