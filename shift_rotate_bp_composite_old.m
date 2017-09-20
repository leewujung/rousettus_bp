function [raw,shift_max,shift_elpctr,shift_elpctr_tilt] = ...
        shift_rotate_bp_composite_old(azq,elq,vq_norm,map_proj,it_shift_th,azel_bnd)
% Shift beampattern according to the center of the best-fitting ellipse and
% compensate for the ellipse rotation

% INPUT
%   azq        azimuth [deg]
%   elq        elevation [deg]
%   vq_norm    normalized bp energy
%   map_proj   map projection to use, 'ortho' or 'eckert4'
%   it_shift_th  threshold for iterative shift ellipse center
%   azel_bnd     bound for az-el shift

% Wu-Jung Lee | leewujung@gmail.com
% 2016 04 19  deal with the "rotate out-of-bound" problem
%             move plotting routine outside of the function to rotate_all_click.m
% 2016 07 21  force E.x0/y0 to be real during iterative rotation
% 2016 07 23  small modification for use with composite clicks

mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);

[xxx,yyy] = mfwdtran(mstruct,[0,0,-90,90],[-180,180,0,0]);
xy_lim = [xxx(1:2), yyy(3:4)];

% [~,vq_norm,azq,elq] = interp_bp(az,el,call_dB,'rbf');  % use the first frequency data for finding ellipse center
% az = az/pi*180;  % convert to [deg]
% el = el/pi*180;
% azq = azq/pi*180;
% elq = elq/pi*180;

% Rotate measurements to use max position as origin
[~,mmidx] = max(vq_norm(:));
origin_max = [elq(mmidx),azq(mmidx)];  % [Lat Lon]
% [el_max,az_max] = rotatem(el,az,origin_max,'forward','degrees');
[elq_max,azq_max] = rotatem(elq,azq,origin_max,'forward','degrees');

% % Determine right/left click
% if azq(mmidx)<0
%     click_side = 0;  % left: click_side = 0
% else
%     click_side = 1;  % right: click_side = 1
% end

% Project lat-lon to map projection distance
% [x,y] = mfwdtran(mstruct,el,az);
[xq,yq] = mfwdtran(mstruct,elq,azq);
% [x_max,y_max] = mfwdtran(mstruct,el_max,az_max);  % map projection
[xq_max,yq_max] = mfwdtran(mstruct,elq_max,azq_max);  % map projection

% raw.x = x;
% raw.y = y;
raw.xq = xq;
raw.yq = yq;
% raw.az = az;
% raw.el = el;
raw.azq = azq;
raw.elq = elq;
raw.vq_norm = vq_norm;
% raw.call_dB = call_dB;

% shift_max.x = x_max;
% shift_max.y = y_max;
shift_max.xq = xq_max;
shift_max.yq = yq_max;
% shift_max.az = az_max;
% shift_max.el = el_max;
shift_max.azq = azq_max;
shift_max.elq = elq_max;
shift_max.vq_norm = vq_norm;
% shift_max.call_dB = call_dB;

E_max = bp_fit_ellipse_azel(shift_max,mstruct);  % fit ellipse
shift_max.E = E_max;


% Iteratively shift until the center of best-fitting ellipse
% is at [0,0] in X-Y plane
n = 0;
E1 = E_max;
% M1.x = x_max;
% M1.y = y_max;
% M1.el = el_max;
% M1.az = az_max;
M1.xq = xq_max;
M1.yq = yq_max;
M1.azq = azq_max;
M1.elq = elq_max;
M1.vq_norm = vq_norm;
% M1.call_dB = call_dB;

% it_shift_th = 0.005;
while sqrt(E1.x0^2+E1.y0^2)>it_shift_th
    [origin1(1),origin1(2)] = ...
        minvtran(mstruct,real(E1.x0),real(E1.y0));  % transform ellipse center x-y pos to lat-lon
    [M1.elq,M1.azq] = rotatem(M1.elq,M1.azq,origin1,'forward','degrees');  % rotation
%     [M1.el,M1.az] = rotatem(M1.el,M1.az,origin1,'forward','degrees');  % rotation
    [M1.xq,M1.yq] = mfwdtran(mstruct,M1.elq,M1.azq);  % map projection
%     [M1.x,M1.y] = mfwdtran(mstruct,M1.el,M1.az);  % map projection
    n = n+1;
    E1 = bp_fit_ellipse_azel(M1,mstruct);
end
fprintf('Rotate original measurements %d times\n',n);
% x_ecen = M1.x;
% y_ecen = M1.y;
xq_ecen = M1.xq;
yq_ecen = M1.yq;

shift_elpctr = M1;

E_ecen = bp_fit_ellipse_azel(shift_elpctr,mstruct);  % fit ellipse
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
% x_ecen_tilt = x_ecen*cos(elp_theta) - y_ecen*sin(elp_theta);
% y_ecen_tilt = x_ecen*sin(elp_theta) + y_ecen*cos(elp_theta);

% Convert from x-y to az-el
% idx_good = x_ecen_tilt>xy_lim(1) & x_ecen_tilt<xy_lim(2) &...
%     y_ecen_tilt>xy_lim(3) & y_ecen_tilt<xy_lim(4);
% el_ecen_tilt = nan(length(x_ecen_tilt),1);
% az_ecen_tilt = nan(length(y_ecen_tilt),1);
% [el_ecen_tilt(idx_good),az_ecen_tilt(idx_good)] = ...
%     minvtran(mstruct,x_ecen_tilt(idx_good),y_ecen_tilt(idx_good));

idx_q_good = xq_ecen_tilt>xy_lim(1) & xq_ecen_tilt<xy_lim(2) &...
    yq_ecen_tilt>xy_lim(3) & yq_ecen_tilt<xy_lim(4);
elq_ecen_tilt = nan(size(xq_ecen_tilt));
azq_ecen_tilt = nan(size(xq_ecen_tilt));
[elq_ecen_tilt(idx_q_good),azq_ecen_tilt(idx_q_good)] = ...
    minvtran(mstruct,xq_ecen_tilt(idx_q_good),yq_ecen_tilt(idx_q_good));

% shift_elpctr_tilt.x = x_ecen_tilt;
% shift_elpctr_tilt.y = y_ecen_tilt;
shift_elpctr_tilt.xq = xq_ecen_tilt;
shift_elpctr_tilt.yq = yq_ecen_tilt;
% shift_elpctr_tilt.az = az_ecen_tilt;
% shift_elpctr_tilt.el = el_ecen_tilt;
shift_elpctr_tilt.azq = azq_ecen_tilt;
shift_elpctr_tilt.elq = elq_ecen_tilt;
shift_elpctr_tilt.vq_norm = vq_norm;
% shift_elpctr_tilt.call_dB = call_dB;

E_ecen_tilt = bp_fit_ellipse_azel(shift_elpctr_tilt,mstruct);  % check tilt compensation results
shift_elpctr_tilt.E = E_ecen_tilt;


% azq, elq shift/tilt out-of-bnd check
% azel_bnd = [300,150];
ofb_max_q = abs(shift_max.azq-raw.azq)>azel_bnd(1) |...
            abs(shift_max.elq-raw.elq)>azel_bnd(2);
ofb_elpctr_q = abs(shift_elpctr.azq-shift_max.azq)>azel_bnd(1) |...
               abs(shift_elpctr.elq-shift_max.elq)>azel_bnd(2) | ofb_max_q;
ofb_tilt_q = abs(shift_elpctr_tilt.azq-shift_elpctr.azq)>azel_bnd(1) |...
             abs(shift_elpctr_tilt.elq-shift_elpctr.elq)>azel_bnd(2) | ofb_elpctr_q;

shift_max.outofbnd_azel_idx_q = ofb_max_q;
shift_elpctr.outofbnd_azel_idx_q = ofb_elpctr_q;
shift_elpctr_tilt.outofbnd_azel_idx_q = ofb_tilt_q;

% az, el shift/tilt check
% ofb_max = abs(shift_max.az-raw.az)>azel_bnd(1) |...
%           abs(shift_max.el-raw.el)>azel_bnd(2);
% ofb_elpctr = abs(shift_elpctr.az-shift_max.az)>azel_bnd(1) |...
%              abs(shift_elpctr.el-shift_max.el)>azel_bnd(2) | ofb_max;
% ofb_tilt = abs(shift_elpctr_tilt.az-shift_elpctr.az)>azel_bnd(1) |...
%            abs(shift_elpctr_tilt.el-shift_elpctr.el)>azel_bnd(2) | ofb_elpctr;

% shift_max.outofbnd_azel_idx = ofb_max;
% shift_elpctr.outofbnd_azel_idx = ofb_elpctr;
% shift_elpctr_tilt.outofbnd_azel_idx = ofb_tilt;



     

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
