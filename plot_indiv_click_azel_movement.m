function fig_azel = plot_indiv_click_azel_movement(Vpass)

% Wu-Jung Lee | leewujung@gmail.com
% 2016/10/25  Do not plot outofbnd_azel_idx az-el locations
%             because this was taken out from shift_rotate_bp
%             since it's reasonable to have az=180 to -180 and el=90 to -90
% 2017/03/08  Change the part on assembling azel for plotting,
%             since the new shift_rotate_bp function would give empty
%             struct if the estimated ellipse center is beyond the
%             boundary of globe (line 82)

fig_azel = figure;
for iP=1:length(Vpass.raw.az)
    azel = [Vpass.raw.az(iP), Vpass.raw.el(iP);...
        Vpass.rot_max.az(iP), Vpass.rot_max.el(iP)];
    if ~isempty(Vpass.rot_elpctr)
        azel = [azel; Vpass.rot_elpctr.az(iP), Vpass.rot_elpctr.el(iP)];
    end
    if ~isempty(Vpass.rot_elpctr_tilt)
        azel = [azel; Vpass.rot_elpctr_tilt.az(iP), Vpass.rot_elpctr_tilt.el(iP)];
    end
    figure(fig_azel)
    plot(azel(:,1),azel(:,2),'.-');  % each line trace the az/el movement of each mic
    hold on
% ====== below commented 2016/10/25 ==============
%     if ismember(iP,find(Vpass.rot_elpctr_tilt.outofbnd_azel_idx))  % highlight the out-of-bound mic
%         plot(azel(:,1),azel(:,2),'ro-','linewidth',2);
%     end
% ================================================
    text(azel(1,1),azel(1,2),num2str(iP));
    axis([-180 180 -90 90])
    grid on
    xlabel('Azimuth (deg)');
    ylabel('Elevation (deg)');
end
