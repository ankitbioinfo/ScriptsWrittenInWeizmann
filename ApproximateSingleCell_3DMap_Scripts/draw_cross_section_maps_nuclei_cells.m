function [] = draw_cross_section_maps_nuclei_cells(centroid_nuc_or_cel,features,grid,visualization,opt,names,f,min_val,max_val,range,cmap,minaxislimit,maxaxislimit)

            
  
    %cross section Z    
    h1 = figure('units','normalized','outerposition',[0 0 1 1]);
    colors_cross_sections = jet(5);
    count=1;
    clear CN_belongs_to_grid
    clear CN_color_belong_to_grid
    for cs =1:5,
        cscount=1;
         for i = 1:grid.nb_x,
            for j = 1:grid.nb_y,
                for k = floor(1+(cs-1)*grid.nb_z/5):ceil(cs*grid.nb_z/5),
                    if features.(names{2}).vals(i,j,k) >=opt.minimum_no_object_in_grid,
                          %[cs,i,j,k]
                          subplot(2,3,cs+1)
                          hold on,
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
                          surface(x,y,z,'FaceColor',  cmap(idx_tmp,:),'EdgeColor',cmap(idx_tmp,:));
                          index_of_CN=find((grid.x_bins(i)<=centroid_nuc_or_cel(:,1))&(centroid_nuc_or_cel(:,1)<grid.x_bins(i+1))& (grid.y_bins(j)<=centroid_nuc_or_cel(:,2))&(centroid_nuc_or_cel(:,2)<grid.y_bins(j+1))& (grid.z_bins(k)<=centroid_nuc_or_cel(:,3))&(centroid_nuc_or_cel(:,3)<grid.z_bins(k+1)));
                          CN_belongs_to_grid{cs+1}{cscount}=index_of_CN;
                          CN_color_belong_to_grid{cs+1}(cscount,:)=cmap(idx_tmp,:);
                          hold off,     
                          
                          subplot(2,3,1)
                          hold on,
                          surface(x,y,z,'FaceColor',  colors_cross_sections(cs,:),'EdgeColor',colors_cross_sections(cs,:))
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
    xlabel('x');ylabel('y');zlabel('z');
    axis image
    title(features.(names{f}).title,'Interpreter','None');
    view([visualization.az visualization.el]);
    if opt.save_figs,
        cs_folder=[opt.save_folder,'/cross_section_features/'];
        if ~exist(cs_folder,'dir')
            mkdir([cs_folder]);
        end
        if opt.matlabfig
		tic
            saveas(h1,[cs_folder,names{f},'_cross_section_Z']);
		
        end
	
        saveas(h1,[cs_folder,names{f},'_cross_section_Z','.png']);
		toc
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
              title(features.(names{f}).title,'Interpreter','None');
        else
              view([0,90]);     
              title(['Cross section ', num2str(cs-1)],'Color', colors_cross_sections(cs-1,:))
        end
    end
    
    if opt.save_figs,
        cs_folder=[opt.save_folder,'/cross_section_features/'];
        if ~exist(cs_folder,'dir')
            mkdir([cs_folder]);
        end
        tic
        if opt.matlabfig
            saveas(h1,[cs_folder,names{f},'_cross_section_Z_patches']);
        end
        saveas(h1,[cs_folder,names{f},'_cross_section_Z_patches','.png']);
        toc
    end
    close all
    
