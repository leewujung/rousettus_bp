function E = bp_fit_ellipse_azel(M,mstruct)
% Find best-fitting ellipse to the -3dB contour
% Have option to plot or not

% INPUT
%   M   struct with x,y,az,el,etc data
%   mstruct    map struct for plotting

% Wu-Jung Lee | leewujung@gmail.com
% 2015/12/25  update to use 'get_main_contour.m' to obtain -3dB contour
% 2016/04/19  move plotting routine outside of the function to rotate_all_click.m
% 2016/07/26  check is c_main empty


% Get -3dB contour
idx_notnan = ~isnan(M.azq);  % index of non-NaN data
azq = M.azq(idx_notnan);
elq = M.elq(idx_notnan);
azq = min(azq(:)):max(azq(:));
elq = min(elq(:)):max(elq(:));
[azq,elq] = meshgrid(azq,elq);
[azq,elq,vq] = griddata(M.azq(idx_notnan),M.elq(idx_notnan),M.vq_norm(idx_notnan),azq,elq);  % griddata to create regular orthogonal azq/elq

[c_main,~] = get_main_contour(vq,unique(azq),unique(elq),-3);  % get main contour at -3dB
if isempty(c_main)
    E = [];
    return;
else
[c3db_xy(:,1),c3db_xy(:,2)] = mfwdtran(mstruct,c_main(:,2),c_main(:,1));  % [az,el] to [x,y]
end

% Ellipse fitting
A = EllipseDirectFit(c3db_xy);  % fit ellipse (direct fit)
E = get_ellipse_param(A);       % get ellipse parameters
E.c3db_xy = c3db_xy;
 

% Old routine for finding -3dB contour: before 2014/12/24
% [C,~] = contour(xq,yq,vq,0:-3:-39,'fill','on');
% Cout = parse_contour_output(C);
% c3db_xy = [];
% for iT=1:length(Cout)  % in case contour break into pieces
%     if Cout(iT).Level == -3
%         c3db_xy = [c3db_xy; Cout(iT).X',Cout(iT).Y'];
%     end
% end
