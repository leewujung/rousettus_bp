function fig_fit = plot_indiv_click_rotate_new(Vpass,cvec,mstruct)
%
% Wu-Jung Lee | leewujung@gmail.com
% 2016/04/19
% 2016/07/21  Add map grid/frame onto the plot results


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
            V = Vpass.rot_elpctr;
            t_txt = sprintf('Tilt %2.2fdeg',V.E.theta/pi*180);
        case 4
            V = Vpass.rot_elpctr_tilt;
            t_txt = 'Tilt compensated';
    end
    
    % valid vq_norm range after rotation
    if iSUB~=1
        V.vq_norm(V.outofbnd_azel_idx_q) = NaN;
    end
    
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

    cvec_new = 0:diff(cvec(1:2)):diff(cvec(1:2))*20;
    %contour(V.xq,V.yq,V.vq_norm,cvec(2:end),'fill','on');
    contourf(V.xq,V.yq,V.vq_norm,cvec_new,'fill','on','linecolor','w');
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
    
    title(t_txt);
    colorbar('ticks',fliplr(cvec(1:2:end)));
    colormap(parula(length(cvec)-1));
    caxis([cvec(end), 0])
end
