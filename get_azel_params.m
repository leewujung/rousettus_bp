function params = get_azel_params(B)
% Obtain az/el-associated values for use in
% `fig_azel_distr_indiv_composite_model_20170923.m`
%
% Wu-Jung Lee | leewujung@gmail.com
% 2017 09 23

if ~isempty(B.rot_elpctr_tilt)
    [ar_tmp,elps_xy_tmp] = get_ellipse_ar(B.rot_elpctr_tilt.E);
    params.ar = ar_tmp;
    params.elps_x = range(elps_xy_tmp(:,1))/2;
    params.elps_y = range(elps_xy_tmp(:,2))/2;
    params.e = B.rot_elpctr_tilt.E.e;
    params.a0 = B.rot_elpctr_tilt.E.a0;
    params.b0 = B.rot_elpctr_tilt.E.b0;
    params.theta = B.rot_elpctr_tilt.E.theta;
else
    params.ar = NaN;
    params.elps_x = NaN;
    params.elps_y  = NaN;
    params.e = NaN;
    params.a0 = NaN;
    params.b0 = NaN;
    params.theta = NaN;
end
