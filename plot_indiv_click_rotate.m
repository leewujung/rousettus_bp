function fig_fit = plot_indiv_click_rotate(Vpass,cvec,mstruct)
%
% Wu-Jung Lee | leewujung@gmail.com
% 2016/04/19
% 2016/07/21  Add map grid/frame onto the plot results
% 2016/10/25  Take out part related to outofbnd_azel_idx_q
% 2017/03/08  Additional check for the potentially empty struct given by 
%             shift_rotate_bp function when the estimated ellipse center
%             is beyond the boundary of globe (line 82)

fig_fit = figure;
set(fig_fit,'position',[140 100 900 480])

for iSUB = 1:4
    switch iSUB
        case 1
            V = Vpass.raw;
            t_txt = 'Raw data';
        case 2
            V = Vpass.rot_max;
            t_txt = 'Best-fitting ellipse';
        case 3
            if ~isempty(Vpass.rot_elpctr)
                V = Vpass.rot_elpctr;
                t_txt = sprintf('Tilt %2.2fdeg',V.E.theta/pi*180);
            else
                continue;
            end
        case 4
            if ~isempty(Vpass.rot_elpctr_tilt)
                V = Vpass.rot_elpctr_tilt;
                t_txt = 'Tilt compensated';
            else
                continue;
            end
    end
    
    % ====== below commented 2016/10/25 ==============
    %     % valid vq_norm range after rotation
    %     if iSUB~=1
    %         V.vq_norm(V.outofbnd_azel_idx_q) = NaN;
    %     end
    % =================================+==============
    
    if ~isempty(V)
        % bound for plotting ellipse
        if iSUB~=1
            E = V.E;
            c3db_xy = E.c3db_xy;
            xmin = min(c3db_xy(:,1))-range(c3db_xy(:,1)*0.5);
            xmax = max(c3db_xy(:,1))+range(c3db_xy(:,1)*0.5);
            ymin = min(c3db_xy(:,2))-range(c3db_xy(:,1)*0.5);
            ymax = max(c3db_xy(:,2))+range(c3db_xy(:,1)*0.5);
        end
        
        % plot contour and ellipse
        subplot(2,2,iSUB);
        
        axesm(mstruct);  % set map axes
        axis off
        
        contour(V.xq,V.yq,V.vq_norm,cvec(2:end),'fill','on');
        if iSUB~=1
            hold on
            fit_df = ezplot(E.eqt,[xmin,xmax,ymin,ymax]);  % plot ellipse
            set(fit_df,'linecolor','b','linewidth',2);
            plot(E.x0,E.y0,'r*');
            text(E.x0,E.y0,sprintf('%2.3f, %2.3f',E.x0,E.y0))
        end
        framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
        gridm('gcolor',200*ones(1,3)/255,'glinestyle','-');
        tightmap
        %axis equal; grid on
        %axis(xy_lim)
        
        title(t_txt);
        %                 colorbar('ticks',fliplr(cvec(1:2:end)));
        colormap(parula(length(cvec)-1));
        caxis([cvec(end), 0])
    end
end
