%% Plot Single Particle Model w/ Electrolyte & Temperature (SPMeT) Results
%   Published December 18, 2016 by Professor Scott Moura
%   Energy, Controls, and Applications Lab (eCAL)
%   University of California, Berkeley
%   http://ecal.berkeley.edu/

close all;
fs = 16;

%% Plot Current, SOC, Voltage
%%% Static Plant Data
figure(1); clf;
set(gcf,'Position',[99,29,623,669],'PaperPositionMode','auto');

% Current
subplot(411);
plot(t,I/p.OneC,'LineWidth',2);
ylabel('Current [C-Rate]','FontSize',fs);
legend({'$$I(t)$$'},'interpreter','latex','Fontsize',fs)
set(gca,'FontSize',fs)
xlim([0, t(end)])
set(gca,'XTickLabel',{''});
set(gca,'Position',[0.12 0.81 0.82 0.18])

% Surface Concentrations
subplot(412);
plot(t,c_ss_n/p.c_s_n_max,'b-',t,c_ss_p/p.c_s_p_max,'r-','LineWidth',2);
hold on;
plot(t,SOC_n,'b--',t,SOC_p,'r--','LineWidth',2);
leg_css = {'$$\theta^-(t)$$';'$$\theta^+(t)$$';'$$\overline{\theta}^-(t)$$';'$$\overline{\theta}^+(t)$$'};
legend(leg_css,'FontSize',fs,'Interpreter','latex','Fontsize',fs)
ylabel('Surf. & Bulk SOC. [-]','FontSize',fs);
set(gca,'FontSize',fs)
xlim([0, t(end)])
set(gca,'XTickLabel',{''});
set(gca,'Position',[0.12 0.57 0.82 0.22])

% Temperature
subplot(413);
plot(t,T1-273.15,'c-',t,T2-273.15,'m-','LineWidth',2);
ylabel('Temperature [C]','FontSize',fs);
legend({'$$T_1(t)$$';'$$T_2(t)$$'},'interpreter','latex','Fontsize',fs)
set(gca,'FontSize',fs)
xlim([0, t(end)])
set(gca,'XTickLabel',{''});
set(gca,'Position',[0.12 0.32 0.82 0.22])

% Voltage
subplot(414);
plot(t,V,'LineWidth',2);
ylabel('Voltage [V]','FontSize',fs);
xlabel('Time [sec]','FontSize',fs);
legend({'$$V(t)$$'},'interpreter','latex','Fontsize',fs)
set(gca,'FontSize',fs)
xlim([0, t(end)])
set(gca,'Position',[0.12 0.08 0.82 0.22])

