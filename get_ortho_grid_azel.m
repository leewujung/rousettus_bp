function [azq,elq,vq] = get_ortho_grid_azel(azq,elq,vq)

idx_notnan = ~isnan(azq);  % index of non-NaN data
azq = azq(idx_notnan);
elq = elq(idx_notnan);
azq = min(azq(:)):max(azq(:));
elq = min(elq(:)):max(elq(:));
[azq,elq] = meshgrid(azq,elq);
[azq,elq,vq] = griddata(azq(idx_notnan),elq(idx_notnan),vq(idx_notnan),azq,elq);  % griddata to create regular orthogonal azq/elq
