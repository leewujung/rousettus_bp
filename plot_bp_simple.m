function h = plot_bp_simple(h,az,el,v,map_proj)
% Utility function to plot beampattern on a globe.
% If az/el are 1D, apply interpolation on v to obtain contourf numbers.
% If az/el are 2D, directly plot using numbers in v.
% Both az and el have unit [deg]
%
% Wu-Jung Lee | leewujung@gmail.com
% 2015 11 04  Test rotating max amplitude point to middle of az-el plane
% 2015 11 20  Use new format in call_dB selection
% 2017 09 21  Morph from plot_bp_on_globe.m

vq_norm_min = -27;
contour_vec = 0:-3:(floor(vq_norm_min/3)-1)*3;
cvec_min_idx = find(contour_vec-vq_norm_min<0,1,'first');

bnd_idx = boundary(az(:),el(:),0);  % outer boundary of all measured points
bnd = [az(bnd_idx),el(bnd_idx)];

if any(size(az)==1)  % az/el/v are vectors
                     % Interpolation
    [azq,elq] = meshgrid(min(az):max(az),min(el):max(el));
    vq = rbfinterp([azq(:)';elq(:)'],rbfcreate([az(:)';el(:)'],v,...
                                               'RBFFunction','multiquadrics'));
    vq = reshape(vq,size(azq));
else  % az/el/v are matrix form, no need for interpolation
    azq = az;
    elq = el;
    vq = v;
end
[in,on] = inpolygon(azq(:),elq(:),bnd(:,1),bnd(:,2));
in_smpl_poly = in|on;
clear in on
vq(~in_smpl_poly) = NaN;
   
% Plot
axes(h)
cla
axesm(map_proj);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
axis off
contourfm(elq,azq,vq,contour_vec(2:cvec_min_idx),...
          'fill','on','linecolor','w');  % don't plot 0 dB contour

colorbar('southoutside','ticks',fliplr(contour_vec(1:cvec_min_idx)));
colormap(parula(cvec_min_idx-1));
caxis([contour_vec(cvec_min_idx) 0]);
tightmap

