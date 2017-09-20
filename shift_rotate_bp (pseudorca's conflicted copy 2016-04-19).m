function [click_side,raw,shift_max,shift_elpctr,shift_elpctr_tilt,fig_fit] = shift_rotate_bp(az,el,call_dB,map_proj,plot_indiv)
% Shift beampattern according to the center of the best-fitting ellipse and
% compensate for the ellipse rotation

% INPUT
%   az   azimuth [radian]
%   el   elevation [radian]
%   call_dB    call energy from each mic [dB]
%   map_proj   map projection to use, 'ortho' or 'eckert4'
%   plot_indiv   if plot individual shift/rotate procedure

% Wu-Jung Lee | leewujung@gmail.com
% 2015 12 27  First stable version
% 2016 04 18  Deal with "rotate out-of-bound" problem

mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);

% Get x-y limit of projected map
[xxx,yyy] = mfwdtran(mstruct,[0,0,-90,90],[-180,180,0,0]);
xy_lim = [xxx(1:2), yyy(3:4)];

[~,vq_norm,azq,elq] = interp_bp(az,el,call_dB,'rbf');  % use the first frequency data for finding ellipse center
az = az/pi*180;  % convert to [deg]
el = el/pi*180;
azq = azq/pi*180;
elq = elq/pi*180;

% Rotate measurements to use max position as origin
[~,mmidx] = max(vq_norm(:));
origin_max = [elq(mmidx),azq(mmidx)];  % [Lat Lon]
[el_max,az_max] = rotatem(el,az,origin_max,'forward','degrees');
[elq_max,azq_max] = rotatem(elq,azq,origin_max,'forward','degrees');

% azel = [az_max-az,el_max-el];
% azelq = [azq_max(:)-azq(:),elq_max(:)-elq(:)];
% azel_outofbnd_idx = ~(abs(azel(:,1))>300 | abs(azel(:,2))>150);  % find out-of-bound rotation
% outofbnd_idxq = ~(abs(azelq(:,1))>300 | abs(azelq(:,2))>150);  % find out-of-bound rotation
% 
% vq_norm_fix = vq_norm;
% vq_norm_fix(~outofbnd_idxq) = NaN;  % interpolated bp without out-of-bnd points

% Determine right/left click
if azq(mmidx)<0
    click_side = 0;  % left: click_side = 0
else
    click_side = 1;  % right: click_side = 1
end

% Project lat-lon to map projection distance
[x,y] = mfwdtran(mstruct,el,az);
[xq,yq] = mfwdtran(mstruct,elq,azq);
[x_max,y_max] = mfwdtran(mstruct,el_max,az_max);  % map projection
[xq_max,yq_max] = mfwdtran(mstruct,elq_max,azq_max);  % map projection
% x_max(~good_idx) = [];  % delete points that are rotated out-of-bnd
% y_max(~good_idx) = [];

% Plot raw data
if plot_indiv==1
    fig_fit = figure;
    set(fig_fit,'position',[140 100 900 480])
    subplot(221)
    contour(xq,yq,vq_norm,-3:-3:-39,'fill','on');
    axis equal; axis(xy_lim); grid on
    title('Raw data');
end

raw.x = x;
raw.y = y;
raw.xq = xq;
raw.yq = yq;
raw.az = az;
raw.el = el;
raw.azq = azq;
raw.elq = elq;
raw.vq_norm = vq_norm;
raw.call_dB = call_dB;

shift_max.x = x_max;
shift_max.y = y_max;
shift_max.xq = xq_max;
shift_max.yq = yq_max;
shift_max.az = az_max;
shift_max.el = el_max;
shift_max.azq = azq_max;
shift_max.elq = elq_max;
shift_max.vq_norm = vq_norm;
shift_max.vq_norm_fix = vq_norm_fix;
shift_max.call_dB = call_dB;

% Plot measurements rotated to max point and fit ellipse
if plot_indiv
    E_max = bp_fit_ellipse_azel(shift_max,mstruct,subplot(222));  % fit ellipse
    title('Best-fitting ellipse');
else
    E_max = bp_fit_ellipse_azel(shift_max,mstruct);  % fit ellipse
end
shift_max.E = E_max;

% Iteratively shift until the center of best-fitting ellipse
% is at [0,0] in X-Y plane
n = 0;
E1 = E_max;
M1.x = x_max;
M1.y = y_max;
M1.el = el_max;
M1.az = az_max;
M1.xq = xq_max;
M1.yq = yq_max;
M1.azq = azq_max;
M1.elq = elq_max;
M1.vq_norm = vq_norm;
M1.call_dB = call_dB;
while sqrt(E1.x0^2+E1.y0^2)>0.005
    [origin1(1),origin1(2)] = minvtran(mstruct,E1.x0,E1.y0);  % transform ellipse center x-y pos to lat-lon
    [M1.elq,M1.azq] = rotatem(M1.elq,M1.azq,origin1,'forward','degrees');  % rotation
    [M1.el,M1.az] = rotatem(M1.el,M1.az,origin1,'forward','degrees');  % rotation
    [M1.xq,M1.yq] = mfwdtran(mstruct,M1.elq,M1.azq);  % map projection
    [M1.x,M1.y] = mfwdtran(mstruct,M1.el,M1.az);  % map projection
    n = n+1;
    E1 = bp_fit_ellipse_azel(M1,mstruct);
end
fprintf('Rotate original measurements %d times\n',n);
x_ecen = M1.x;
y_ecen = M1.y;
xq_ecen = M1.xq;
yq_ecen = M1.yq;

shift_elpctr = M1;

if plot_indiv
    E_ecen = bp_fit_ellipse_azel(shift_elpctr,mstruct,subplot(223));  % fit ellipse
    title(sprintf('Tilt %2.2fdeg',E_ecen.theta/pi*180));
else
    E_ecen = bp_fit_ellipse_azel(shift_elpctr,mstruct);  % fit ellipse
end
shift_elpctr.E = E_ecen;


% Compensate for ellipse tilt
if E_ecen.theta<0 && E_ecen.coef.a<E_ecen.coef.c && E_ecen.a0>E_ecen.b0  % horizontal ellipse
    elp_theta = -(pi/2+E_ecen.theta);
elseif E_ecen.theta>0 && E_ecen.coef.a<E_ecen.coef.c && E_ecen.a0>E_ecen.b0  % horizontal ellipse
    elp_theta = pi/2-E_ecen.theta;
else
    elp_theta = -E_ecen.theta;
end

xq_ecen_tilt = xq_ecen*cos(elp_theta) - yq_ecen*sin(elp_theta);
yq_ecen_tilt = xq_ecen*sin(elp_theta) + yq_ecen*cos(elp_theta);
x_ecen_tilt = x_ecen*cos(elp_theta) - y_ecen*sin(elp_theta);
y_ecen_tilt = x_ecen*sin(elp_theta) + y_ecen*cos(elp_theta);

% Convert from x-y to az-el
idx_good = x_ecen_tilt>xy_lim(1) & x_ecen_tilt<xy_lim(2) & y_ecen_tilt>xy_lim(3) & y_ecen_tilt<xy_lim(4);
el_ecen_tilt = nan(length(x_ecen_tilt),1);
az_ecen_tilt = nan(length(y_ecen_tilt),1);
[el_ecen_tilt(idx_good),az_ecen_tilt(idx_good)] = minvtran(mstruct,x_ecen_tilt(idx_good),y_ecen_tilt(idx_good));

idx_q_good = xq_ecen_tilt>xy_lim(1) & xq_ecen_tilt<xy_lim(2) & yq_ecen_tilt>xy_lim(3) & yq_ecen_tilt<xy_lim(4);
elq_ecen_tilt = nan(size(xq_ecen_tilt));
azq_ecen_tilt = nan(size(xq_ecen_tilt));
[elq_ecen_tilt(idx_q_good),azq_ecen_tilt(idx_q_good)] = minvtran(mstruct,xq_ecen_tilt(idx_q_good),yq_ecen_tilt(idx_q_good));

shift_elpctr_tilt.x = x_ecen_tilt;
shift_elpctr_tilt.y = y_ecen_tilt;
shift_elpctr_tilt.xq = xq_ecen_tilt;
shift_elpctr_tilt.yq = yq_ecen_tilt;
shift_elpctr_tilt.az = az_ecen_tilt;
shift_elpctr_tilt.el = el_ecen_tilt;
shift_elpctr_tilt.azq = azq_ecen_tilt;
shift_elpctr_tilt.elq = elq_ecen_tilt;
shift_elpctr_tilt.vq_norm = vq_norm;
shift_elpctr_tilt.call_dB = call_dB;

if plot_indiv
    E_ecen_tilt = bp_fit_ellipse_azel(shift_elpctr_tilt,mstruct,subplot(224));  % check tilt compensation results
    title('Tilt compensated')
else
    E_ecen_tilt = bp_fit_ellipse_azel(shift_elpctr_tilt,mstruct);  % check tilt compensation results
end
shift_elpctr_tilt.E = E_ecen_tilt;


if plot_indiv==0
    fig_fit = [];
end



% If center of the rotated ellipse is not close to [0,0]
% if sqrt(shift_elpctr_tilt.E.x0^2+shift_elpctr_tilt.E.y0^2)>0.005
%     n = 0;
%     E2 = E_ecen_tilt;
%     M2 = shift_elpctr_tilt;
% %     fig_test = figure;
% %     haxes = axes;
%     while sqrt(E2.x0^2+E2.y0^2)>0.005
%         M2_tmp = M2;
%         [origin1(1),origin1(2)] = minvtran(mstruct,E2.x0,E2.y0);  % transform ellipse center x-y pos to lat-lon
%         [M2_tmp.elq,M2_tmp.azq] = rotatem(M2_tmp.elq,M2_tmp.azq,origin1,'forward','degrees');  % rotation
%         [M2_tmp.el,M2_tmp.az] = rotatem(M2_tmp.el,M2_tmp.az,origin1,'forward','degrees');  % rotation
%         [M2_tmp.xq,M2_tmp.yq] = mfwdtran(mstruct,M2_tmp.elq,M2_tmp.azq);  % map projection
%         [M2_tmp.x,M2_tmp.y] = mfwdtran(mstruct,M2_tmp.el,M2_tmp.az);  % map projection
%         n = n+1;
%         
%         E2_tmp = bp_fit_ellipse_azel(M2,mstruct);
% %         E2_tmp = bp_fit_ellipse_azel(M2_tmp,mstruct,haxes);
%         if sqrt(E2_tmp.x0^2+E2_tmp.y0^2)>=sqrt(E2.x0^2+E2.y0^2)
%             break
%         else
%             E2 = E2_tmp;
%             M2 = M2_tmp;
%         end
%     end
%     fprintf('Rotate rotated results %d times\n',n);
%     shift_elpctr_tilt_new = M2;
%     
%     if plot_indiv
%         E_ecen_tilt_new = bp_fit_ellipse_azel(shift_elpctr_tilt_new,mstruct,subplot(235));  % Fit ellipse
%         title('Shift to [0,0] again');
%     else
%         E_ecen_tilt_new = bp_fit_ellipse_azel(shift_elpctr_tilt_new,mstruct);  % Fit ellipse
%     end
%     shift_elpctr_tilt_new.E = E_ecen_tilt_new;
%     
%     shift_elpctr_tilt = shift_elpctr_tilt_new;  % Update shift-and-tilt results
% end
