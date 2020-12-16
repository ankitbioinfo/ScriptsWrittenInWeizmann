function [] = draw_maps_nuclei_cells_DT_and_PT(visualization,opt,type,b)
%DRAW_MAPS Draws 3D maps from grid features

    if exist([opt.path{b},'Alignment_matrix.dat'],'file')==0
        registration=eye(3);
    else
 	M = load([opt.path{b},'Alignment_matrix.dat']);
 	registration=M(:,1:3);
    end
 	mainVec=registration(:,3); 
        clear M
 	 

disp('draw_maps_nuclei_cells_DT_and_PT')
for f = 2:12,
   
    if type == 1
        res_saved = load([opt.path_PT{b},'nuclei_grid.mat']);
        features = res_saved.nuclei_grid;
        grid = res_saved.grid;
    elseif type == 2
        res_saved = load([opt.path_PT{b},'cells_grid.mat']);
        features = res_saved.cells_grid;
        grid = res_saved.grid;
    end
    names = fieldnames(features);
    names{f}
    min_val = visualization.(names{f}).min;
    max_val = visualization.(names{f}).max;
    range = round(min_val*features.(names{f}).scaling_factor):round(max_val*features.(names{f}).scaling_factor);
    cmap = jet(length(range));
    
    h1 = figure(1);
    
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
    
    if type == 1
        res_saved = load([opt.path_DT{b},'nuclei_grid.mat']);
        features = res_saved.nuclei_grid;
        grid = res_saved.grid;
    elseif type == 2
        res_saved = load([opt.path_DT{b},'cells_grid.mat']);
        features = res_saved.cells_grid;
        grid = res_saved.grid;
    end
    
    names = fieldnames(features);
    
    hold on
    for i = 1:grid.nb_x,
        for j = 1:grid.nb_y,
            for k = 1:grid.nb_z,
                if features.(names{2}).vals(i,j,k) >= 10,
                    r = (opt.delta_x/2)*(0.8*((features.(names{f}).vals(i,j,k) - min_val)/...
                        (max_val - min_val))-0.2*((features.(names{f}).vals(i,j,k) - max_val)/(max_val - min_val)));
                    
                    [xs,ys,zs] = sphere(50);
                    
                    x0 = features.centers(i,j,k,1)+opt.DT_shift_x;
                    y0 = features.centers(i,j,k,2)+opt.DT_shift_y;
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

% Log Volume and Surface Area
disp('log volume and surface area')

for f = 3:4
    
    if type == 1
        res_saved = load([opt.path_PT{b},'nuclei_grid.mat']);
        features = res_saved.nuclei_grid;
        grid = res_saved.grid;
    elseif type == 2
        res_saved = load([opt.path_PT{b},'cells_grid.mat']);
        features = res_saved.cells_grid;
        grid = res_saved.grid;
    end
    
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
    
    if type == 1
        res_saved = load([opt.path_DT{b},'nuclei_grid.mat']);
        features = res_saved.nuclei_grid;
        grid = res_saved.grid;
    elseif type == 2
        res_saved = load([opt.path_DT{b},'cells_grid.mat']);
        features = res_saved.cells_grid;
        grid = res_saved.grid;
    end
    
    hold on
    
    for i = 1:grid.nb_x,
        for j = 1:grid.nb_y,
            for k = 1:grid.nb_z,
                if features.(names{2}).vals(i,j,k) >=10,
                    
                    r = (opt.delta_x/2)*(0.8*((features.(names{f}).scaling_factor_log*log(features.(names{f}).vals(i,j,k)) - min_val)/(max_val - min_val))-0.2*((features.(names{f}).scaling_factor_log*log(features.(names{f}).vals(i,j,k)) - max_val)/(max_val - min_val)));
                    [xs,ys,zs] = sphere(50);
                    
                    x0 = features.centers(i,j,k,1)+opt.DT_shift_x;
                    y0 = features.centers(i,j,k,2)+opt.DT_shift_y;
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
    %xlim([7000 9000])
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

% Orientations
disp('Orientation PCs')
for f = (length(names)-3):(length(names)-1)
    
    if type == 1
        res_saved = load([opt.path_PT{b},'nuclei_grid.mat']);
        features = res_saved.nuclei_grid;
        grid = res_saved.grid;
    elseif type == 2
        res_saved = load([opt.path_PT{b},'cells_grid.mat']);
        features = res_saved.cells_grid;
        grid = res_saved.grid;
    end
    
    names = fieldnames(features);
    
    amp = 60;
    h1 = figure(1);

    THmin = 1;
    THmax = 1;
    PHImin = 1;
    PHImax = 1;

    for i = 1:grid.nb_x,
        for j = 1:grid.nb_y,
            for k = 1:grid.nb_z,
                if features.(names{2}).vals(i,j,k) >=10,
                    hold on,
%                     TH = acos(features.(names{f}).orientations(i,j,k,3));
%                     PHI = atan2(features.(names{f}).orientations(i,j,k,2),features.(names{f}).orientations(i,j,k,1));
%                     THmin  = min(THmin, TH);
%                     PHImin  = min(THmin, TH);
%                     THmax  = max(THmax, TH);
%                     PHImax  = max(THmax, TH);
                      xvec= [features.centers(i,j,k,1);features.centers(i,j,k,1)+amp*features.(names{f}).orientations(i,j,k,1)];
                      yvec= [features.centers(i,j,k,2);features.centers(i,j,k,2)+amp*features.(names{f}).orientations(i,j,k,2)];
                      zvec= [features.centers(i,j,k,3);features.centers(i,j,k,3)+amp*features.(names{f}).orientations(i,j,k,3)];
                      paulVector=[diff(xvec),diff(yvec),diff(zvec)];                          
                      theta=angleCompute(paulVector,mainVec);
                      %mycolor=[(1/2*sin(TH*2) + 1/2),(1/2*cos(PHI) + 1/2),(1/2*sin(PHI) + 1/2)];
                      line(xvec,yvec,zvec,'LineWidth',15*features.(names{f}).sph_var(i,j,k),'Color',pickcolor(theta))
                      hold off,
                end
            end
        end
    end
    
    if type == 1
        res_saved = load([opt.path_DT{b},'nuclei_grid.mat']);
        features = res_saved.nuclei_grid;
        grid = res_saved.grid;
    elseif type == 2
        res_saved = load([opt.path_DT{b},'cells_grid.mat']);
        features = res_saved.cells_grid;
        grid = res_saved.grid;
    end
    
    names = fieldnames(features);
    
    hold on
    
    for i = 1:grid.nb_x,
        for j = 1:grid.nb_y,
            for k = 1:grid.nb_z,
                if features.(names{2}).vals(i,j,k) >=10,

                    hold on,
%                     TH = acos(features.(names{f}).orientations(i,j,k,3));
%                     PHI = atan2(features.(names{f}).orientations(i,j,k,2),features.(names{f}).orientations(i,j,k,1));
%                     THmin  = min(THmin, TH);
%                     PHImin  = min(THmin, TH);
%                     THmax  = max(THmax, TH);
%                     PHImax  = max(THmax, TH);
                      xvec= [features.centers(i,j,k,1);features.centers(i,j,k,1)+amp*features.(names{f}).orientations(i,j,k,1)];
                      yvec= [features.centers(i,j,k,2);features.centers(i,j,k,2)+amp*features.(names{f}).orientations(i,j,k,2)];
                      zvec= [features.centers(i,j,k,3);features.centers(i,j,k,3)+amp*features.(names{f}).orientations(i,j,k,3)];
                      paulVector=[diff(xvec),diff(yvec),diff(zvec)];                          
                      theta=angleCompute(paulVector,mainVec);
                      %mycolor=[(1/2*sin(TH*2) + 1/2),(1/2*cos(PHI) + 1/2),(1/2*sin(PHI) + 1/2)];
                      line(xvec,yvec,zvec,'LineWidth',15*features.(names{f}).sph_var(i,j,k),'Color',pickcolor(theta))
                      hold off,
                end
            end
        end
    end

    xlabel('x');
    ylabel('y');
    zlabel('z');

    axis image
    title(names{f})

    view([visualization.az visualization.el]);

    if opt.save_figs,
        if ~exist([opt.save_folder],'dir')
            mkdir([opt.save_folder]);
        end
        if opt.matlabfig
            saveas(h1,[opt.save_folder,names{f},'_orientation']);
        end
        saveas(h1,[opt.save_folder,names{f},'_orientation','.png']);
    end
    close all
end

end


 %angleCompute=@(u,v) atan2(norm(cross(u,v)),dot(u,v));
function value=angleCompute(u,v) 
         value=atan2(norm(cross(u,v)),dot(u,v));
end



function chooseColor=pickcolor(value)
        colormap default;
        if value<=pi/2
           value=value;
        elseif (value<=pi)&(value>pi/2)
           value=pi-value;
        elseif (value<=1.5*pi)&(value>pi)
           value=value-pi;
        elseif (value<=2*pi)&(value>1.5*pi)
           value=2*pi-value;
        end
        
        map=colormap;
        theta=linspace(0,pi/2,65);
        for i=2:length(theta)
            if value<theta(i)
                chooseColor=map(i-1,:);
                break 
            end
        end
 end
    
