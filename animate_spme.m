%% Animate Single Particle Model w/ Electrolyte (SPMe) Results
%   Published November 5, 2016 by Professor Scott Moura
%   Energy, Controls, and Applications Lab (eCAL)
%   University of California, Berkeley
%   http://ecal.berkeley.edu/

close all;
fs = 16;

%     vidObj = VideoWriter('img/SPMe_UDDS.avi');
%     vidObj.FrameRate = 10;
%     vidObj.Quality = 100;
%     open(vidObj);
    
figure(2);
set(gcf,'Position',[461    13   803   693],'PaperPositionMode','auto');

cemaxmax = max(max(c_e));
ceminmin = min(min(c_e));

for k = 1:NT

    clf;

    % Anode Concentration
    subplot(3,3,1)
    plot(r_vec,c_n(k,:)/p.c_s_n_max,'LineWidth',2);
    ylim([0,1])
    ylabel('Solid Conc., \theta^-(t)','FontSize',fs);
    xlabel('Radial coordinate, r [node no.]','FontSize',fs)
    title('\bf ANODE','fontsize',fs)
    set(gca,'FontSize',fs);
    set(gca,'Position',[0.1 0.70 0.28 0.27])

    % Cathode Concentration
    subplot(3,3,3)
    plot(r_vec(end:-1:1),c_p(k,:)/p.c_s_p_max,'LineWidth',2);
    ylim([0,1])
    ylabel('Solid Conc., \theta^+(t)','FontSize',fs);
    xlabel('Radial coordinate, r [node no.]','FontSize',fs)
    title('\bf CATHODE','fontsize',fs)
    set(gca,'FontSize',fs);
    set(gca,'Position',[0.67 0.7 0.28 0.27])
    XTickLabel = get(gca,'XTickLabel');
    set(gca,'XTickLabel',XTickLabel(end:-1:1,:));

    % Electrolyte Concentration
    subplot(3,3,[4 5 6])
    plot(x_vec_spme,c_e(:,k)/1e3,'b-','LineWidth',2)
    xlim([x_vec_spme(1), x_vec_spme(end)]);
    ylim([ceminmin, cemaxmax]/1e3);
    ylabel('Elec. Conc., c_e(x,t) [kmol/m^3]','FontSize',fs);
    xlabel('Space across cell, x [node no.]','FontSize',fs)
    set(gca,'FontSize',fs);
    set(gca,'XTickLabel','');
    set(gca,'Position',[0.1 0.37 0.85 0.25])

    % Voltage and Current
    subplot(3,3,[7 8 9]);
%     plot(t(1:k),I(1:k),'g-',t(k),I(k),'go','LineWidth',2);
    hold on;
    plot(t(1:k),V(1:k),'b-',t(k),V(k),'bo','LineWidth',2,'MarkerSize',10);
    xlim([0, t(end)])
    ylim([min(V), max(V)]);
    ylabel('Voltage [V]','FontSize',fs)
    xlabel('Time [sec]','FontSize',fs);
    set(gca,'FontSize',fs);
    set(gca,'Position',[0.1 0.08 0.85 0.24])

%         F = getframe(gcf);
%         writeVideo(vidObj,F);

    pause(0.05);

end

% close(vidObj);