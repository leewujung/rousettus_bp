function [c_level,c_level_nan] = get_main_contour(vq,azq,elq,level)

% INPUT
%   vq     interpolated beampattern values
%   azq    azimuth [degree]
%   elq    azimuth [degree]
%   level  level of beampattern value wanted
% OUTPUT
%   c_level        [az,el] contour at the level of specified bp value
%   c_level_wnan   same as c_level but with NaN inserted at measurement boundary

% Wu-Jung Lee | leewujung@gmail.com
% 2015/12/24
% 2016/07/26  Check if max_idx is empty --> when all contours are outside
%             of globe

% Find contour surrounding the area with measurement
vq_new = zeros(size(vq));
vq_new(~isnan(vq)) = 1;
c_edge = contourc(azq,elq,vq_new,1);
c_edge = parse_contour_output(c_edge);
c_edge = [c_edge.X; c_edge.Y]';

% Find contour at the specified level
im_cc = bwconncomp(vq>=level,4);
im_area = regionprops(im_cc,'Area');
[~,max_idx] = max([im_area(:).Area]);
if isempty(max_idx)
    c_level = [];
    c_level_nan = [];
    return;
end
im_main_idx = im_cc.PixelIdxList{max_idx};
im_main = zeros(size(vq));
im_main(im_main_idx) = 1;
im_main = imfill(im_main);
c_level = contourc(azq,elq,im_main,1);
c_level = parse_contour_output(c_level);
c_level = [c_level.X; c_level.Y]';

% fig_chk = figure;
% plot(c_edge(:,1),c_edge(:,2));
% hold on
% plot(c_level(:,1),c_level(:,2));

% Delete overlapping portion
D = pdist2(c_level,c_edge);
[idx_level_zero,~] = ind2sub(size(D),find(D==0));
c_level_nan = c_level;
c_level(idx_level_zero,:) = [];
c_level_nan(idx_level_zero,:) = NaN;

% figure(fig_chk)
% plot(c_level_nan(:,1),c_level_nan(:,2),'linewidth',1.5);
% disp('checked')

