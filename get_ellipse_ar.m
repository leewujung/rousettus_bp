function [ar,elps_xy] = get_ellipse_ar(E)

c3db_xy = E.c3db_xy;
xmin = min(c3db_xy(:,1))-range(c3db_xy(:,1)*0.5);
xmax = max(c3db_xy(:,1))+range(c3db_xy(:,1)*0.5);
ymin = min(c3db_xy(:,2))-range(c3db_xy(:,1)*0.5);
ymax = max(c3db_xy(:,2))+range(c3db_xy(:,1)*0.5);

fit_df = ezplot(E.eqt,[xmin,xmax,ymin,ymax]);  % plot ellipse
tmp = get(fit_df,'contourMatrix');

% best-fittign ellipse in eckert 4 projection
elps_xy(:,1) = tmp(1,2:end)';
elps_xy(:,2) = tmp(2,2:end)';

% aspect ratio
ar = range(elps_xy(:,2))/range(elps_xy(:,1));

close