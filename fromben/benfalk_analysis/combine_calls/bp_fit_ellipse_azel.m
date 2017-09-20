function E = bp_fit_ellipse_azel(M,mstruct,varargin)
% Find best-fitting ellipse to the -3dB contour
% Have option to plot or not

% INPUT
%   M   struct with x,y,az,el,etc data
%   mstruct    map struct for plotting
%   varargin   haxes, if exist then plot

% Wu-Jung Lee | leewujung@gmail.com
% 2015/12/25  update to use 'get_main_contour.m' to obtain -3dB contour

% Parse input
if nargin>2
    plot_opt = 1;
    haxes = varargin{1};
else
    plot_opt = 0;
end

% Get -3dB contour
idx_notnan = ~isnan(M.azq);  % index of non-NaN data
azq = M.azq(idx_notnan);
elq = M.elq(idx_notnan);
azq = min(azq(:)):max(azq(:));
elq = min(elq(:)):max(elq(:));
[azq,elq] = meshgrid(azq,elq);
[azq,elq,vq] = griddata(M.azq(idx_notnan),M.elq(idx_notnan),M.vq_norm(idx_notnan),azq,elq);  % griddata to create regular orthogonal azq/elq

[c_main,~] = get_main_contour(vq,unique(azq),unique(elq),-3);  % get main contour at -3dB
[c3db_xy(:,1),c3db_xy(:,2)] = mfwdtran(mstruct,c_main(:,2),c_main(:,1));  % [az,el] to [x,y]

% Ellipse fitting
A = EllipseDirectFit(c3db_xy);  % fit ellipse (direct fit)
E = get_ellipse_param(A);       % get ellipse parameters
xmin = min(c3db_xy(:,1))-range(c3db_xy(:,1)*0.5);
xmax = max(c3db_xy(:,1))+range(c3db_xy(:,1)*0.5);
ymin = min(c3db_xy(:,2))-range(c3db_xy(:,1)*0.5);
ymax = max(c3db_xy(:,2))+range(c3db_xy(:,1)*0.5);

% Plot
if plot_opt==1
    % find x-y plot limit
    [xxx,yyy] = mfwdtran(mstruct,[0,0,-90,90],[-180,180,0,0]);
    xy_lim = [xxx(1:2), yyy(3:4)];
    
    % contour plot with best-fitting ellipse overlaid
    axes(haxes);
    contour(M.xq,M.yq,M.vq_norm,-39:3:-3,'fill','on');
    hold on
    fit_df = ezplot(E.eqt,[xmin,xmax,ymin,ymax]);  % plot ellipse
    set(fit_df,'linecolor','b','linewidth',2);
    plot(E.x0,E.y0,'r*');
    text(E.x0,E.y0,sprintf('%2.3f, %2.3f',E.x0,E.y0))
    title('')
    hold off
    axis equal; grid on
    axis(xy_lim);
end


% Old routine for finding -3dB contour: before 2014/12/24
% [C,~] = contour(xq,yq,vq,0:-3:-39,'fill','on');
% Cout = parse_contour_output(C);
% c3db_xy = [];
% for iT=1:length(Cout)  % in case contour break into pieces
%     if Cout(iT).Level == -3
%         c3db_xy = [c3db_xy; Cout(iT).X',Cout(iT).Y'];
%     end
% end
