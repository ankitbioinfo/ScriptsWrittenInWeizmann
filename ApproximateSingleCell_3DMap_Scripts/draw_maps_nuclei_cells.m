function [] = draw_maps_nuclei_cells(cel_nuc_data,features,grid,visualization,opt,gp,cellornuclei)
%DRAW_MAPS Draws 3D maps from grid features
%
  names = fieldnames(features); 
  centroid_nuc_or_cel=cel_nuc_data(:,5:7);
%   1  'centers'
%   2  'object_number'
%   3  'volume'
%   4  'surface_area'
%   5  'sphericity'
%   6  'Allometric'
%   7  'occupation'
%   8  'density'
%   9  'density_delaunay'
%   10  'ratio_PC3_PC2'
%   11  'ratio_PC2_PC1'
%   12  'ratio_PC3_PC1'
%   13  'PC1'
%   14  'PC2'
%   15  'PC3'
  
%  M = load([opt.path{gp},'Alignment_matrix.dat']);
%  vec=M(:,1:3);val=M(:,4:6);

 
disp('draw_maps_nuclei_cells')
if opt.plot_3D_map_14_features==1
for f = 10:length(names),
    names{f}
    min_val = visualization.(names{f}).min;
    max_val = visualization.(names{f}).max;
    range = round(min_val*features.(names{f}).scaling_factor):round(max_val*features.(names{f}).scaling_factor);
    cmap = jet(length(range));
    
    h1 = figure(1);
    
    count=1;
    clear minaxislimit
    clear maxaxislimit 
    for i = 1:grid.nb_x,
        for j = 1:grid.nb_y,
            for k = 1:grid.nb_z,
                if features.(names{2}).vals(i,j,k) >= opt.minimum_no_object_in_grid,
                    r = (opt.delta_x/2)*(0.8*((features.(names{f}).vals(i,j,k) - min_val)/...
                        (max_val - min_val))-0.2*((features.(names{f}).vals(i,j,k) - max_val)/(max_val - min_val)));
                    
                    [xs,ys,zs] = sphere(50);
                    
                    x0 = features.centers(i,j,k,1);
                    y0 = features.centers(i,j,k,2);
                    z0 = features.centers(i,j,k,3);
                    
                    x = xs*r + x0;
                    y = ys*r + y0;
                    z = zs*r + z0;
                    
                    minaxislimit(count,:)=[x0-r,y0-r,z0-r];
                    maxaxislimit(count,:)=[x0+r,y0+r,z0+r];
                    count=count+1;
                    
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
    
    draw_cross_section_maps_nuclei_cells(centroid_nuc_or_cel,features,grid,visualization,opt,names,f,min_val,max_val,range,cmap,minaxislimit,maxaxislimit)
    
    
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
end




if opt.plot_3D_map_log_features==1
%Log Volume and Surface Area
disp('log volume and surface area')
for f = 3:4
    min_val = features.(names{f}).scaling_factor_log*log(visualization.(names{f}).min);
    max_val = features.(names{f}).scaling_factor_log*log(visualization.(names{f}).max);
    range = round(min_val):round(max_val);
    cmap = jet(length(range));
    
    h1 = figure(1);
    
    for i = 1:grid.nb_x,
        for j = 1:grid.nb_y,
            for k = 1:grid.nb_z,
                if features.(names{2}).vals(i,j,k) >=opt.minimum_no_object_in_grid,
                    
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


%Calculate average orientations of bone 
%The angle deviation in next orientation map used colorVec  
    for f = (length(names)-2):(length(names))
         names{f};
         length(names);
         count=1;
         for i = 1:grid.nb_x,
            for j = 1:grid.nb_y,
                   for k = 1:grid.nb_z,
                    if features.(names{2}).vals(i,j,k) >=opt.minimum_no_object_in_grid,
                          xvec(count,1)= [features.(names{f}).orientations(i,j,k,1)];
                          yvec(count,1)= [features.(names{f}).orientations(i,j,k,2)];
                          zvec(count,1)= [features.(names{f}).orientations(i,j,k,3)];       
                          count=count+1;
                    end
                end
            end
         end
         vals=[xvec yvec zvec];
         vals = [vals; -vals];
         [coeff,score,latent] = pca(vals);
          indtemp2 = find(dot(vals(:,:)',repmat(coeff(:,1),[1,size(vals(:,:),1)]))>0);
          R = sqrt((sum(vals(indtemp2,1))).^2+(sum(vals(indtemp2,2))).^2+(sum(vals(indtemp2,3))).^2);
          colorVec(:,f) = [sum(vals(indtemp2,1))/R, sum(vals(indtemp2,2))/R,sum(vals(indtemp2,3))/R];   
          clear xvec
          clear yvec
          clear zvec
    end

if opt.plot_PC_orientations==1
disp('Orientation PCs')
for f = (length(names)-2):(length(names)) 
    mainVec=colorVec(:,f); 
    amp = 60;
    h1 = figure(1);
    count=1;
    clear data
    clear CN_belongs_to_grid
    clear CN_color_belong_to_grid
    for i = 1:grid.nb_x,
        for j = 1:grid.nb_y,
            for k = 1:grid.nb_z,
                if features.(names{2}).vals(i,j,k) >=opt.minimum_no_object_in_grid,
                    hold on,
                      xvec= [features.centers(i,j,k,1);features.centers(i,j,k,1)+amp*features.(names{f}).orientations(i,j,k,1)];
                      yvec= [features.centers(i,j,k,2);features.centers(i,j,k,2)+amp*features.(names{f}).orientations(i,j,k,2)];
                      zvec= [features.centers(i,j,k,3);features.centers(i,j,k,3)+amp*features.(names{f}).orientations(i,j,k,3)];
                      paulVector=[diff(xvec),diff(yvec),diff(zvec)];                          
                      theta=angleCompute(paulVector,mainVec);
                      data(count,:)=[features.centers(i,j,k,1),features.centers(i,j,k,2),features.centers(i,j,k,3),theta*180/pi,convertAngle(theta)*180/pi];
                      %mycolor=[(1/2*sin(TH*2) + 1/2),(1/2*cos(PHI) + 1/2),(1/2*sin(PHI) + 1/2)];
                      if opt.minimum_no_object_in_grid==1
                          lw=0.5;
                      else
                          lw=15*features.(names{f}).sph_var(i,j,k);
                      end
                      line(xvec,yvec,zvec,'LineWidth',lw,'Color',pickcolor(theta));
                      CN_belongs_to_grid{count}=find((grid.x_bins(i)<=centroid_nuc_or_cel(:,1))&(centroid_nuc_or_cel(:,1)<grid.x_bins(i+1))& (grid.y_bins(j)<=centroid_nuc_or_cel(:,2))&(centroid_nuc_or_cel(:,2)<grid.y_bins(j+1))& (grid.z_bins(k)<=centroid_nuc_or_cel(:,3))&(centroid_nuc_or_cel(:,3)<grid.z_bins(k+1)));
                      CN_color_belong_to_grid(count,:)=pickcolor(theta);
                                    
                      hold off,        
                      count=count+1;
                end
            end
        end
    end
    
    orientationData{f}=data;
    
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
        saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_orientation']);
        end
        saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_orientation','.png']);
    end
    
    close all    
    h1 = figure(1);
    [~,ia]=unique(CN_color_belong_to_grid,'rows');
    uniqueColor=CN_color_belong_to_grid(ia,:);    
    for j=1:size(uniqueColor,1)
        combined_cells_of_grids=[];
        for i=1:size(CN_color_belong_to_grid,1)
             if sum(CN_color_belong_to_grid(i,:)==uniqueColor(j,:))==3
                 combined_cells_of_grids=[combined_cells_of_grids;CN_belongs_to_grid{i}];
             end
        end
        index=combined_cells_of_grids;
        c=centroid_nuc_or_cel(index,:); shp = alphaShape(c);
        plot(shp,'facecolor',uniqueColor(j,:),'edgecolor','none','facealpha',0.8); 
        hold on 
    end   
    xlabel('x');ylabel('y');zlabel('z');
    view([visualization.az visualization.el]);     
    if opt.save_figs,
        if ~exist([opt.save_folder],'dir')
            mkdir([opt.save_folder]);
        end
        if opt.matlabfig
        saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_orientation_patches']);
        end
        saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_orientation_patches','.png']);
    end
    close all 
      
    
end
  

% Orientation Range 
map=flipud(colormap(jet(7)));
theta=[0,5,10,15,25,40,60,90];
h1=figure();
set(gcf, 'PaperSize', [2 5]);
set(gcf, 'PaperPosition', [0 0 2 5]);
tname={'<5','5-10','10-15','15-25','25-40','40-60','60-90'};

for i=2:length(theta)
plot(0.1,-i*0.7,'marker','o','markersize',20, 'markerfacecolor',map(i-1,:),'color',map(i-1,:))
    hold on;
    text(0.3,-i*0.7,tname{i-1},'fontsize',18);
end

axis([0,1,-6.2,-0.8])
title('Orientation colormap')

set(findobj(gcf,'type','axes'),'visible','off');

    if opt.save_figs,
        if ~exist([opt.save_folder],'dir')
            mkdir([opt.save_folder]);
        end
        if opt.matlabfig
        saveas(h1,[opt.save_folder,cellornuclei,'_PCs_orientation_range']);
        end
        saveas(h1,[opt.save_folder,cellornuclei,'_PCs_orientation_range','.png']);
    end

close all 
end
    



if opt.plot_PC_cross_section_orientationZ==1
% Orientations cross sections Z
disp('Orientation Cross sections along Z')
for f = (length(names)-2):(length(names))  
    mainVec=colorVec(:,f);             
    amp = 60;
    h1 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,3,1)    
    count=1;
    clear minaxislimit
    clear maxaxislimit 
    clear CN_belongs_to_grid
    clear CN_color_belong_to_grid
    for i = 1:grid.nb_x,
        for j = 1:grid.nb_y,
            for k = 1:grid.nb_z,
                if features.(names{2}).vals(i,j,k) >=opt.minimum_no_object_in_grid,
                     hold on,
                       xvec= [features.centers(i,j,k,1);features.centers(i,j,k,1)+amp*features.(names{f}).orientations(i,j,k,1)];
                       yvec= [features.centers(i,j,k,2);features.centers(i,j,k,2)+amp*features.(names{f}).orientations(i,j,k,2)];  
                       zvec= [features.centers(i,j,k,3);features.centers(i,j,k,3)+amp*features.(names{f}).orientations(i,j,k,3)];                              
                          
                      minaxislimit(count,:)=[min(xvec),min(yvec),min(zvec)];
                      maxaxislimit(count,:)=[max(xvec),max(yvec),max(zvec)];
                      paulVector=[diff(xvec),diff(yvec),diff(zvec)];                          
                      theta=angleCompute(paulVector,mainVec);
                      count=count+1;    
                      %mycolor=[(1/2*sin(TH*2) + 1/2),(1/2*cos(PHI) + 1/2),(1/2*sin(PHI) + 1/2)];
                      if opt.minimum_no_object_in_grid==1
                          lw=0.5;
                      else
                          lw=15*features.(names{f}).sph_var(i,j,k);
                      end
                      line(xvec,yvec,zvec,'LineWidth',lw,'Color',pickcolor(theta))
                      hold off,
                end
            end
        end
    end
    
    xlabel('x');ylabel('y');zlabel('z');
    axis image
    title(names{f})
    view([visualization.az visualization.el]);
    colors_cross_sections = jet(5);
    count=1;
    for cs =1:5,
        cscount=1;
        %paul cut around X axis 
        %for i = floor(1+(cs-1)*grid.nb_x/5):ceil(cs*grid.nb_x/5),
         for i = 1:grid.nb_x,
            for j = 1:grid.nb_y,
                %for k = 1:grid.nb_z,  cut around Z axis 
                for k = floor(1+(cs-1)*grid.nb_z/5):ceil(cs*grid.nb_z/5),
                    if features.(names{2}).vals(i,j,k) >=opt.minimum_no_object_in_grid,
                          subplot(2,3,cs+1)
                          hold on,
                           xvec= [features.centers(i,j,k,1);features.centers(i,j,k,1)+amp*features.(names{f}).orientations(i,j,k,1)];
                           yvec= [features.centers(i,j,k,2);features.centers(i,j,k,2)+amp*features.(names{f}).orientations(i,j,k,2)];  
                           zvec= [features.centers(i,j,k,3);features.centers(i,j,k,3)+amp*features.(names{f}).orientations(i,j,k,3)];                              
                                
                          paulVector=[diff(xvec),diff(yvec),diff(zvec)];                          
                          theta=angleCompute(paulVector,mainVec); 
                          if opt.minimum_no_object_in_grid==1
                            lw=0.5;
                          else
                            lw=15*features.(names{f}).sph_var(i,j,k);
                          end
                          line(xvec,yvec,zvec,'LineWidth',lw,'Color',pickcolor(theta))
                          index_of_CN=find((grid.x_bins(i)<=centroid_nuc_or_cel(:,1))&(centroid_nuc_or_cel(:,1)<grid.x_bins(i+1))& (grid.y_bins(j)<=centroid_nuc_or_cel(:,2))&(centroid_nuc_or_cel(:,2)<grid.y_bins(j+1))& (grid.z_bins(k)<=centroid_nuc_or_cel(:,3))&(centroid_nuc_or_cel(:,3)<grid.z_bins(k+1)));
                          CN_belongs_to_grid{cs+1}{cscount}=index_of_CN;
                          CN_color_belong_to_grid{cs+1}(cscount,:)=pickcolor(theta);
                         
                          axis image
                          hold off,     
                          
                          subplot(2,3,1)
                          hold on,
                          line(xvec,yvec,zvec,'LineWidth',lw,'Color',colors_cross_sections(cs,:))
                          CN_belongs_to_grid{1}{count}=index_of_CN;
                          CN_color_belong_to_grid{1}(count,:)=colors_cross_sections(cs,:);
                          count=count+1;
                          cscount=cscount+1;
                          hold off,
                    end
                end
            end
        end
        subplot(2,3,cs+1)
        xlabel('x');
        ylabel('y');
        zlabel('z');
        axis([min(minaxislimit(:,1)), max(maxaxislimit(:,1)), min(minaxislimit(:,2)),max(maxaxislimit(:,2)), min(minaxislimit(:,3)),max(maxaxislimit(:,3))]);
        title(['Cross section ', num2str(cs)],'Color', colors_cross_sections(cs,:))
        view([0,90])       
    end
    subplot(2,3,1)
    view([visualization.az visualization.el]);
    if opt.save_figs,
        if ~exist([opt.save_folder],'dir')
            mkdir([opt.save_folder]);
        end
        if opt.matlabfig
            saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_cross_section_orientationZ']);
        end
        saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_cross_section_orientationZ','.png']);
    end
    close all
    
    %cross section patche mode 
    h1 = figure('units','normalized','outerposition',[0 0 1 1]);    
    for cs=1:6 
        subplot(2,3,cs)        
        [~,ia]=unique(CN_color_belong_to_grid{cs},'rows');
        uniqueColor=CN_color_belong_to_grid{cs}(ia,:);    
        for j=1:size(uniqueColor,1)
            combined_cells_of_grids=[];
            for i=1:size(CN_color_belong_to_grid{cs},1)
                 if sum(CN_color_belong_to_grid{cs}(i,:)==uniqueColor(j,:))==3
                     combined_cells_of_grids=[combined_cells_of_grids;CN_belongs_to_grid{cs}{i}];
                 end
            end
            index=combined_cells_of_grids;
            c=centroid_nuc_or_cel(index,:);   shp = alphaShape(c);
            plot(shp,'facecolor',uniqueColor(j,:),'edgecolor','none','facealpha',0.8); 
            hold on 
        end
    
        xlabel('x');ylabel('y');zlabel('z');
        axis image
        axis([min(minaxislimit(:,1)), max(maxaxislimit(:,1)), min(minaxislimit(:,2)),max(maxaxislimit(:,2)), min(minaxislimit(:,3)),max(maxaxislimit(:,3))]);
        if cs==1
              view([visualization.az visualization.el]);
              title(names{f})
        else
              view([0,90]);     
              title(['Cross section ', num2str(cs-1)],'Color', colors_cross_sections(cs-1,:))
        end
    end
    
    if opt.save_figs,
        if ~exist([opt.save_folder],'dir')
            mkdir([opt.save_folder]);
        end
        if opt.matlabfig
            saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_cross_section_orientationZ_patches']);
        end
        saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_cross_section_orientationZ_patches','.png']);
    end
    close all
end     
end  
    

 


if opt.plot_PC_cross_section_orientationX==1
% Orientations cross sections X
disp('Orientation Cross sections along X')
for f = (length(names)-2):(length(names))  
    mainVec=colorVec(:,f);             
    amp = 60;
    h1 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,3,1)    
    count=1;
    clear minaxislimit
    clear maxaxislimit 
    clear CN_belongs_to_grid
    clear CN_color_belong_to_grid
    for i = 1:grid.nb_x,
        for j = 1:grid.nb_y,
            for k = 1:grid.nb_z,
                if features.(names{2}).vals(i,j,k) >=opt.minimum_no_object_in_grid,
                     hold on,
                       xvec= [features.centers(i,j,k,1);features.centers(i,j,k,1)+amp*features.(names{f}).orientations(i,j,k,1)];
                       yvec= [features.centers(i,j,k,2);features.centers(i,j,k,2)+amp*features.(names{f}).orientations(i,j,k,2)];  
                       zvec= [features.centers(i,j,k,3);features.centers(i,j,k,3)+amp*features.(names{f}).orientations(i,j,k,3)];                              
                          
                      minaxislimit(count,:)=[min(xvec),min(yvec),min(zvec)];
                      maxaxislimit(count,:)=[max(xvec),max(yvec),max(zvec)];
                      paulVector=[diff(xvec),diff(yvec),diff(zvec)];                          
                      theta=angleCompute(paulVector,mainVec);
                      count=count+1;    
                      %mycolor=[(1/2*sin(TH*2) + 1/2),(1/2*cos(PHI) + 1/2),(1/2*sin(PHI) + 1/2)];
                      if opt.minimum_no_object_in_grid==1
                          lw=0.5;
                      else
                          lw=15*features.(names{f}).sph_var(i,j,k);
                      end
                      line(xvec,yvec,zvec,'LineWidth',lw,'Color',pickcolor(theta))
                      hold off,
                end
            end
        end
    end
    
    xlabel('x');ylabel('y');zlabel('z');
    axis image
    title(names{f})
    view([135 10]);
    colors_cross_sections = jet(5);
    count=1;
    for cs =1:5,
        cscount=1;
        %paul cut around X axis 
         for i = floor(1+(cs-1)*grid.nb_x/5):ceil(cs*grid.nb_x/5),
         %for i = 1:grid.nb_x,
            for j = 1:grid.nb_y,
                for k = 1:grid.nb_z,   
                %for k = floor(1+(cs-1)*grid.nb_z/5):ceil(cs*grid.nb_z/5),
                    if features.(names{2}).vals(i,j,k) >=opt.minimum_no_object_in_grid,
                          subplot(2,3,cs+1)
                          hold on,
                           xvec= [features.centers(i,j,k,1);features.centers(i,j,k,1)+amp*features.(names{f}).orientations(i,j,k,1)];
                           yvec= [features.centers(i,j,k,2);features.centers(i,j,k,2)+amp*features.(names{f}).orientations(i,j,k,2)];  
                           zvec= [features.centers(i,j,k,3);features.centers(i,j,k,3)+amp*features.(names{f}).orientations(i,j,k,3)];                              
                                
                          paulVector=[diff(xvec),diff(yvec),diff(zvec)];                          
                          theta=angleCompute(paulVector,mainVec);  
                          if opt.minimum_no_object_in_grid==1
                            lw=0.5;
                          else
                            lw=15*features.(names{f}).sph_var(i,j,k);
                          end
                          line(xvec,yvec,zvec,'LineWidth',lw,'Color',pickcolor(theta))
                          index_of_CN=find((grid.x_bins(i)<=centroid_nuc_or_cel(:,1))&(centroid_nuc_or_cel(:,1)<grid.x_bins(i+1))& (grid.y_bins(j)<=centroid_nuc_or_cel(:,2))&(centroid_nuc_or_cel(:,2)<grid.y_bins(j+1))& (grid.z_bins(k)<=centroid_nuc_or_cel(:,3))&(centroid_nuc_or_cel(:,3)<grid.z_bins(k+1)));
                          CN_belongs_to_grid{cs+1}{cscount}=index_of_CN;
                          CN_color_belong_to_grid{cs+1}(cscount,:)=pickcolor(theta);
                         
                          axis image
                          hold off,     
                          
                          subplot(2,3,1)
                          hold on,
                          line(xvec,yvec,zvec,'LineWidth',lw,'Color',colors_cross_sections(cs,:))
                          CN_belongs_to_grid{1}{count}=index_of_CN;
                          CN_color_belong_to_grid{1}(count,:)=colors_cross_sections(cs,:);
                          count=count+1;
                          cscount=cscount+1;
                          hold off,
                    end
                end
            end
        end
        subplot(2,3,cs+1)
        xlabel('x');
        ylabel('y');
        zlabel('z');
        axis([min(minaxislimit(:,1)), max(maxaxislimit(:,1)), min(minaxislimit(:,2)),max(maxaxislimit(:,2)), min(minaxislimit(:,3)),max(maxaxislimit(:,3))]);
        title(['Cross section ', num2str(cs)],'Color', colors_cross_sections(cs,:))
        view([90,0])       
    end
    subplot(2,3,1)
    view([135 10]);
    if opt.save_figs,
        if ~exist([opt.save_folder],'dir')
            mkdir([opt.save_folder]);
        end
        if opt.matlabfig
            saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_cross_section_orientationX']);
        end
        saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_cross_section_orientationX','.png']);
    end
    close all
    
    %cross section patche mode 
    h1 = figure('units','normalized','outerposition',[0 0 1 1]);    
    for cs=1:6 
        subplot(2,3,cs)        
        [~,ia]=unique(CN_color_belong_to_grid{cs},'rows');
        uniqueColor=CN_color_belong_to_grid{cs}(ia,:);    
        for j=1:size(uniqueColor,1)
            combined_cells_of_grids=[];
            for i=1:size(CN_color_belong_to_grid{cs},1)
                 if sum(CN_color_belong_to_grid{cs}(i,:)==uniqueColor(j,:))==3
                     combined_cells_of_grids=[combined_cells_of_grids;CN_belongs_to_grid{cs}{i}];
                 end
            end
            index=combined_cells_of_grids;
            c=centroid_nuc_or_cel(index,:);   shp = alphaShape(c);
            plot(shp,'facecolor',uniqueColor(j,:),'edgecolor','none','facealpha',0.8); 
            hold on 
        end
    
        xlabel('x');ylabel('y');zlabel('z');
        axis image
        axis([min(minaxislimit(:,1)), max(maxaxislimit(:,1)), min(minaxislimit(:,2)),max(maxaxislimit(:,2)), min(minaxislimit(:,3)),max(maxaxislimit(:,3))]);
        if cs==1
              view([135 10]);
              title(names{f})
        else
              view([90,0]);     
              title(['Cross section ', num2str(cs-1)],'Color', colors_cross_sections(cs-1,:))
        end
    end
    
    if opt.save_figs,
        if ~exist([opt.save_folder],'dir')
            mkdir([opt.save_folder]);
        end
        if opt.matlabfig
            saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_cross_section_orientationX_patches']);
        end
        saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_cross_section_orientationX_patches','.png']);
    end
    close all
end 
end
    
if opt.plot_PC_cross_section_orientationY==1
% Orientations cross sections Y
disp('Orientation Cross sections along Y')
for f = (length(names)-2):(length(names))  
    mainVec=colorVec(:,f);             
    amp = 60;
    h1 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,3,1)    
    count=1;
    clear minaxislimit
    clear maxaxislimit 
    clear CN_belongs_to_grid
    clear CN_color_belong_to_grid
    for i = 1:grid.nb_x,
        for j = 1:grid.nb_y,
            for k = 1:grid.nb_z,
                if features.(names{2}).vals(i,j,k) >=opt.minimum_no_object_in_grid,
                     hold on,
                       xvec= [features.centers(i,j,k,1);features.centers(i,j,k,1)+amp*features.(names{f}).orientations(i,j,k,1)];
                       yvec= [features.centers(i,j,k,2);features.centers(i,j,k,2)+amp*features.(names{f}).orientations(i,j,k,2)];  
                       zvec= [features.centers(i,j,k,3);features.centers(i,j,k,3)+amp*features.(names{f}).orientations(i,j,k,3)];                              
                          
                      minaxislimit(count,:)=[min(xvec),min(yvec),min(zvec)];
                      maxaxislimit(count,:)=[max(xvec),max(yvec),max(zvec)];
                      paulVector=[diff(xvec),diff(yvec),diff(zvec)];                          
                      theta=angleCompute(paulVector,mainVec);
                      count=count+1;    
                      %mycolor=[(1/2*sin(TH*2) + 1/2),(1/2*cos(PHI) + 1/2),(1/2*sin(PHI) + 1/2)];
                      if opt.minimum_no_object_in_grid==1
                          lw=0.5;
                      else
                          lw=15*features.(names{f}).sph_var(i,j,k);
                      end
                      line(xvec,yvec,zvec,'LineWidth',lw,'Color',pickcolor(theta))
                      hold off,
                end
            end
        end
    end
    
    xlabel('x');ylabel('y');zlabel('z');
    axis image
    title(names{f})
    view([75 15]);
    colors_cross_sections = jet(5);
    count=1;
    for cs =1:5,
        cscount=1;
         for i = 1:grid.nb_x,
            for j = floor(1+(cs-1)*grid.nb_y/5):ceil(cs*grid.nb_y/5), 
                for k = 1:grid.nb_z,  
                    if features.(names{2}).vals(i,j,k) >=opt.minimum_no_object_in_grid,
                          subplot(2,3,cs+1)
                          hold on,
                           xvec= [features.centers(i,j,k,1);features.centers(i,j,k,1)+amp*features.(names{f}).orientations(i,j,k,1)];
                           yvec= [features.centers(i,j,k,2);features.centers(i,j,k,2)+amp*features.(names{f}).orientations(i,j,k,2)];  
                           zvec= [features.centers(i,j,k,3);features.centers(i,j,k,3)+amp*features.(names{f}).orientations(i,j,k,3)];                              
                                
                            paulVector=[diff(xvec),diff(yvec),diff(zvec)];                          
                            theta=angleCompute(paulVector,mainVec);
                            if opt.minimum_no_object_in_grid==1
                                lw=0.5;
                            else
                                lw=15*features.(names{f}).sph_var(i,j,k);
                            end
                          line(xvec,yvec,zvec,'LineWidth',lw,'Color',pickcolor(theta))
                          index_of_CN=find((grid.x_bins(i)<=centroid_nuc_or_cel(:,1))&(centroid_nuc_or_cel(:,1)<grid.x_bins(i+1))& (grid.y_bins(j)<=centroid_nuc_or_cel(:,2))&(centroid_nuc_or_cel(:,2)<grid.y_bins(j+1))& (grid.z_bins(k)<=centroid_nuc_or_cel(:,3))&(centroid_nuc_or_cel(:,3)<grid.z_bins(k+1)));
                          CN_belongs_to_grid{cs+1}{cscount}=index_of_CN;
                          CN_color_belong_to_grid{cs+1}(cscount,:)=pickcolor(theta);
                         
                          axis image
                          hold off,     
                          
                          subplot(2,3,1)
                          hold on,
                          line(xvec,yvec,zvec,'LineWidth',lw,'Color',colors_cross_sections(cs,:))
                          CN_belongs_to_grid{1}{count}=index_of_CN;
                          CN_color_belong_to_grid{1}(count,:)=colors_cross_sections(cs,:);
                          count=count+1;
                          cscount=cscount+1;
                          hold off,
                    end
                end
            end
        end
        subplot(2,3,cs+1)
        xlabel('x');
        ylabel('y');
        zlabel('z');
        axis([min(minaxislimit(:,1)), max(maxaxislimit(:,1)), min(minaxislimit(:,2)),max(maxaxislimit(:,2)), min(minaxislimit(:,3)),max(maxaxislimit(:,3))]);
        title(['Cross section ', num2str(cs)],'Color', colors_cross_sections(cs,:))
        view([0,0])       
    end
    subplot(2,3,1)
    view([75 15]);
    if opt.save_figs,
        if ~exist([opt.save_folder],'dir')
            mkdir([opt.save_folder]);
        end
        if opt.matlabfig
            saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_cross_section_orientationY']);
        end
        saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_cross_section_orientationY','.png']);
    end
    close all
    
    %cross section patche mode 
    h1 = figure('units','normalized','outerposition',[0 0 1 1]);    
    for cs=1:6 
        subplot(2,3,cs)        
        [~,ia]=unique(CN_color_belong_to_grid{cs},'rows');
        uniqueColor=CN_color_belong_to_grid{cs}(ia,:);    
        for j=1:size(uniqueColor,1)
            combined_cells_of_grids=[];
            for i=1:size(CN_color_belong_to_grid{cs},1)
                 if sum(CN_color_belong_to_grid{cs}(i,:)==uniqueColor(j,:))==3
                     combined_cells_of_grids=[combined_cells_of_grids;CN_belongs_to_grid{cs}{i}];
                 end
            end
            index=combined_cells_of_grids;
            c=centroid_nuc_or_cel(index,:);   shp = alphaShape(c);
            plot(shp,'facecolor',uniqueColor(j,:),'edgecolor','none','facealpha',0.8); 
            hold on 
        end
    
        xlabel('x');ylabel('y');zlabel('z');
        axis image
        axis([min(minaxislimit(:,1)), max(maxaxislimit(:,1)), min(minaxislimit(:,2)),max(maxaxislimit(:,2)), min(minaxislimit(:,3)),max(maxaxislimit(:,3))]);
        if cs==1
              view([75 15]);
              title(names{f})
        else
              view([0,0]);     
              title(['Cross section ', num2str(cs-1)],'Color', colors_cross_sections(cs-1,:))
        end
    end
    
    if opt.save_figs,
        if ~exist([opt.save_folder],'dir')
            mkdir([opt.save_folder]);
        end
        if opt.matlabfig
            saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_cross_section_orientationY_patches']);
        end
        saveas(h1,[opt.save_folder,cellornuclei,'_',names{f},'_cross_section_orientationY_patches','.png']);
    end
    close all
end
end


end
    
 

%angleCompute=@(u,v) atan2(norm(cross(u,v)),dot(u,v));
function value=angleCompute(u,v) 
         value=atan2(norm(cross(u,v)),dot(u,v));
         value=value*180/pi;
end


function chooseColor=pickcolor(value)
        colormap default;
        if value<=90
           value=value;
        elseif (value<=180)&(value>90)
           value=180-value;
        elseif (value<=270)&(value>180)
           value=value-180;
        elseif (value<=360)&(value>270)
           value=360-value;
        end
        
        map=flipud(colormap(jet(7)));
        theta=[0,5,10,15,25,40,60,90];
        for i=2:length(theta)
            if value<=theta(i)
                chooseColor=map(i-1,:);
                break 
            end
        end
end
    
 
function value=convertAngle(value)
        if value<=pi/2
           value=value;
        elseif (value<=pi)&(value>pi/2)
           value=pi-value;
        elseif (value<=1.5*pi)&(value>pi)
           value=value-pi;
        elseif (value<=2*pi)&(value>1.5*pi)
           value=2*pi-value;
        end
end



