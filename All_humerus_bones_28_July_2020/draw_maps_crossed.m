function [] = draw_maps_crossed(features,grid,visualization,opt)
%DRAW_MAPS Draws 3D maps from grid features
%
names = fieldnames(features);

for f = 2:8,
    
    min_val = visualization.(names{f}).min;
    max_val = visualization.(names{f}).max;
    range = round(min_val*features.(names{f}).scaling_factor):round(max_val*features.(names{f}).scaling_factor);
    cmap = jet(length(range));
    
    h1 = figure(1);
    names{f}
    for i = 1:grid.nb_x,
        for j = 1:grid.nb_y,
            for k = 1:grid.nb_z,
                if features.(names{2}).vals(i,j,k) >= 10,
                    r = (opt.delta_x/2)*(0.8*((features.(names{f}).vals(i,j,k) - min_val)/...
                        (max_val - min_val))-0.2*((features.(names{f}).vals(i,j,k) - max_val)/(max_val - min_val)));
                    
                    [xs,ys,zs] = sphere(50);
                    
                    x0 = features.centers(i,j,k,1);
                    y0 = features.centers(i,j,k,2);
                    z0 = features.centers(i,j,k,3);
                    
                    x = xs*r + x0;
                    y = ys*r + y0;
                    z = zs*r + z0;
                    
                    idx_tmp = find(range == round(features.(names{f}).scaling_factor*features.(names{f}).vals(i,j,k)));
                    
                    surface(x,y,z,'FaceColor',  cmap(idx_tmp,:),'EdgeColor',cmap(idx_tmp,:))
                    hold on
                end
            end
        end
    end
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis image
    title(features.(names{f}).title,'Interpreter','None');
    view([visualization.az visualization.el]);
    
    if opt.save_figs,
        if ~exist([opt.save_folder],'dir')
            mkdir([opt.save_folder]);
        end
          if opt.matlabfig
        saveas(h1,[opt.save_folder,names{f}]);
        end
        saveas(h1,[opt.save_folder,names{f},'.png']);
    end
    close all
    
    
    
    h2 = figure(2);
    axis off;
    for i = 1:7,
        
        val_temp = (max_val - min_val)*(i-1)/6 + min_val;
        r = (opt.delta_x/2)*(0.8*((val_temp - min_val)/(max_val - min_val))-0.2*((val_temp - max_val)/(max_val - min_val)));
        
        [xs,ys,zs] = sphere(50);
        
        x0 = (i-1)*opt.delta_x;
        y0 = 0;
        z0 = 0;
        
        x = xs*r + x0;
        y = ys*r + y0;
        z = zs*r + z0;
        
        idx_tmp = find(range == round(features.(names{f}).scaling_factor*val_temp));
        
        surface(x,y,z,'FaceColor',  cmap(idx_tmp,:),'EdgeColor',cmap(idx_tmp,:));
        hold on
        text(x0-opt.delta_x/3, -opt.delta_y, num2str(val_temp,3),'FontSize',14);
        hold on,
        
    end
    title([features.(names{f}).title,' Range'],'Interpreter','None')
    axis image
    
    if opt.save_figs,
        if ~exist([opt.save_folder],'dir')
            mkdir([opt.save_folder]);
        end
          if opt.matlabfig
        saveas(h2,[opt.save_folder,names{f},'_range']);
        end
        saveas(h2,[opt.save_folder,names{f},'_range','.png']);
    end
    close all
end

% log of centroid shift
for f = 5,
    min_val = features.(names{f}).scaling_factor_log*log(visualization.(names{f}).min);
    max_val = features.(names{f}).scaling_factor_log*log(visualization.(names{f}).max);
    range = round(min_val):round(max_val);
    cmap = jet(length(range));
    
    h1 = figure(1);
    
    for i = 1:grid.nb_x,
        for j = 1:grid.nb_y,
            for k = 1:grid.nb_z,
                if features.(names{2}).vals(i,j,k) >=10,
                    
                    r = (opt.delta_x/2)*(0.8*((features.(names{f}).scaling_factor_log*log(features.(names{f}).vals(i,j,k)) - min_val)/(max_val - min_val))-0.2*((features.(names{f}).scaling_factor_log*log(features.(names{f}).vals(i,j,k)) - max_val)/(max_val - min_val)));
                    [xs,ys,zs] = sphere(50);
                    
                    x0 = features.centers(i,j,k,1);
                    y0 = features.centers(i,j,k,2);
                    z0 = features.centers(i,j,k,3);
                    
                    x = xs*r + x0;
                    y = ys*r + y0;
                    z = zs*r + z0;
                    
                    idx_tmp = find(range == round(features.(names{f}).scaling_factor_log*log(features.(names{f}).vals(i,j,k))));
                    
                    surface(x,y,z,'FaceColor',  cmap(idx_tmp,:),'EdgeColor',cmap(idx_tmp,:))
                    hold on
                end
                
            end
        end
    end
    
    xlabel('x');
    ylabel('y');
    zlabel('z');
    xlim([7000 9000])
    axis image
    title(['Log ',features.(names{f}).title])
    view([visualization.az visualization.el]);
    
    if opt.save_figs,
        if ~exist([opt.save_folder],'dir')
             mkdir([opt.save_folder]);
        end
          if opt.matlabfig
        saveas(h1,[opt.save_folder,names{f},'_log']);
            end
        saveas(h1,[opt.save_folder,names{f},'_log','.png']);
    end
    close all
    
    
    h2 = figure(2);
    axis off;
    for i = 1:7,
        
        val_temp = (max_val - min_val)*(i-1)/6 + min_val;
        r = (opt.delta_x/2)*(0.8*((val_temp - min_val)/(max_val - min_val))-0.2*((val_temp - max_val)/(max_val - min_val)));
        
        [xs,ys,zs] = sphere(50);
        
        x0 = (i-1)*opt.delta_x;
        y0 = 0;
        z0 = 0;
        
        x = xs*r + x0;
        y = ys*r + y0;
        z = zs*r + z0;
        
        idx_tmp = find(range == round(val_temp));
        
        surface(x,y,z,'FaceColor',  cmap(idx_tmp,:),'EdgeColor',cmap(idx_tmp,:));
        hold on
        text(x0-opt.delta_x/3, -opt.delta_y, num2str(exp(val_temp/features.(names{f}).scaling_factor_log),3),'FontSize',14);
        hold on,
        
    end
    title([features.(names{f}).title,' Range'])
    axis image
    
    if opt.save_figs,
        if ~exist([opt.save_folder],'dir')
              mkdir([opt.save_folder]);
        end
          if opt.matlabfig
        saveas(h2,[opt.save_folder,names{f},'_log_range']);
        end
        saveas(h2,[opt.save_folder,names{f},'_log_range','.png']);
    end
    close all

end

end

