function [] = draw_spatial_profiles_new(features,bonetype,grid,opt)
%DRAW_SPATIAL_PROFILES Draws spatial profiles
%

%     'centers'
%     'object_number'
%     'volume'
%     'surface_area'
%     'sphericity'
%     'occupation'
%     'density'
%     'density_delaunay'
%     'PC1'
%     'PC2'
%     'PC3'
%     'ratio_PC1_PC2'

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

names = fieldnames(features);
units = {'','','\mum^3','\mum^2','','','','[per \mum^3]','[per 75 \mum^3]','','','','','',''};

for f = [2:length(names)]

    h = figure('units','normalized','outerposition',[0 0 0.5 1]);
    count=1;
    clear centroid
    clear data 
    for i = 1:grid.nb_x,
            for j = 1:grid.nb_y,
                for k = 1:grid.nb_z,
                    if features.(names{2}).vals(i,j,k) >= opt.minimum_no_object_in_grid,    
                        data(count,:)=features.(names{f}).vals(i,j,k);
                        centroid(count,:)=[features.centers(i,j,k,1),features.centers(i,j,k,2),features.centers(i,j,k,3)];
                        count=count+1;                        
                    end
            end
        end
    end
    
        if (bonetype==1)
            centroidZ=-centroid(:,3);
        else
            centroidZ=centroid(:,3);
        end
        
               
        if strmatch(names{f},'density_delaunay')
            data=data*75^3;   %delaunay density multipled with 75 micron cube that means how many 
                   % cell/nuclei present in 75^3 cube 
        end
        
       
        
        
    minl=min(centroidZ); maxl=max(centroidZ);
    Yinterval=linspace(minl,maxl,opt.grid_bin);  
    %Xinterval=zeros(1,length(Yinterval));
    binsize=Yinterval(2)-Yinterval(1);
    horizontalCount=cell(1,length(Yinterval));
    
    onelessYinterval=Yinterval(1:end-1);
    
       
    clear Xinterval 
    
    for k=1:length(onelessYinterval)
        horizontalCount{k}=[];
        for tt=1:length(centroidZ)
            if (centroidZ(tt) > Yinterval(k)) & (centroidZ(tt) < Yinterval(k+1))
                horizontalCount{k}=[horizontalCount{k},tt];
            end
        end
        Xinterval(k)=mean(data(horizontalCount{k}));
    end
    clear horizontalCount 
    
    plot(centroidZ,data,'b.','markersize',5)
    hold on 
   
    plot(onelessYinterval+binsize/2,Xinterval,'ro-');
    set(gca,'fontsize',10);
    xlabel('Bone long axis RZ-HZ')
    ylabel(['Average ',features.(names{f}).title, ' ',units{f}])

   
    
    
%     valtemp = squeeze(features.(names{f}).vals(:,2,2));
%     indtemp = find(valtemp);
%     plot(features.(names{f}).vals(indtemp,2,2),grid.x_bins(indtemp),'LineWidth',5,'Color','k');
  %  xlim([min(features.(names{f}).vals(indtemp,2,2)) max(features.(names{f}).vals(indtemp,2,2))]);
    title(['Grid average: ', features.(names{f}).title ],'Interpreter','None');
    if opt.save_figs,
        if ~exist([opt.save_folder],'dir')
           mkdir([opt.save_folder]);
        end
        saveas(h,[opt.save_folder,names{f},'_spatial_profile']);
        saveas(h,[opt.save_folder,names{f},'_spatial_profile','.png']);
    end
    close all
    
end

% drawing log spatial profiles for volume and surface area or centroid 
% shift
if isfield(features,'centroid_shift')
    fi = 5;
else
    fi = (3:4); 
end

for f = fi
     h = figure('units','normalized','outerposition',[0 0 0.5 1]);
     
      count=1;
    clear centroid
    clear data 
    for i = 1:grid.nb_x,
            for j = 1:grid.nb_y,
                for k = 1:grid.nb_z,
                    if features.(names{2}).vals(i,j,k) >= opt.minimum_no_object_in_grid,    
                        data(count,:)=features.(names{f}).vals(i,j,k);
                        centroid(count,:)=[features.centers(i,j,k,1),features.centers(i,j,k,2),features.centers(i,j,k,3)];
                        count=count+1;                        
                    end
            end
        end
    end
    
    
        if (bonetype==1)
            centroidZ=-centroid(:,3);
        else
            centroidZ=centroid(:,3);
        end
    minl=min(centroidZ); maxl=max(centroidZ);
    
      
    
    Yinterval=linspace(minl,maxl,opt.grid_bin);  
    %Xinterval=zeros(1,length(Yinterval));
    binsize=Yinterval(2)-Yinterval(1);
    horizontalCount=cell(1,length(Yinterval));
    onelessYinterval=Yinterval(1:end-1);

    clear Xinterval
    for k=1:length(onelessYinterval)
        horizontalCount{k}=[];
        for tt=1:length(centroidZ)
            if (centroidZ(tt) > Yinterval(k)) & (centroidZ(tt) < Yinterval(k+1))
                horizontalCount{k}=[horizontalCount{k},tt];
            end
        end
        Xinterval(k)=mean(data(horizontalCount{k}));
    end
    clear horizontalCount 
    
    semilogy(centroidZ,data,'b.','markersize',5)
    hold on 
    semilogy(onelessYinterval+binsize/2,Xinterval,'ro-');
    set(gca,'fontsize',10);
    xlabel('Bone long axis RZ-HZ')
    ylabel(['average ','log(',features.(names{f}).title,')'])

         
     
%     valtemp = squeeze(features.(names{f}).vals(:,2,2));
%     indtemp = find(valtemp);
%     semilogx(features.(names{f}).vals(indtemp,2,2),grid.x_bins(indtemp),'LineWidth',5,'Color','k');
%     xlim([min(log(features.(names{f}).vals(indtemp,2,2))) max(log(features.(names{f}).vals(indtemp,2,2)))]);
    title(['Grid average: ','log(',features.(names{f}).title,')'],'Interpreter','None');
    if opt.save_figs,
        if ~exist([opt.save_folder],'dir')
            mkdir([opt.save_folder]);
        end
        saveas(h,[opt.save_folder,names{f},'_spatial_profile_log']);
        saveas(h,[opt.save_folder,names{f},'_spatial_profile_log','.png']);
    end
    close all
end

% superimposing volume, surface area and sphericity
% non log and log
% if ~isfield(features,'centroid_shift')
%         h = figure('units','normalized','outerposition',[0 0 0.5 1]);
%     for f = 3:5
%         valtemp = squeeze(features.(names{f}).vals(:,2,2));
%     indtemp = find(valtemp);
%     hold on,
%     plot(features.(names{f}).vals(indtemp,2,2)/max(features.(names{f}).vals(indtemp,2,2)),grid.x_bins(indtemp),'LineWidth',5);
%     hold off,
%     end
%     legend({'Vol','Surf','Spher'});
%     title('Volume, surface area, sphericity','Interpreter','None');
%     if opt.save_figs,
%         if ~exist([opt.save_folder],'dir')
%                 mkdir([opt.save_folder]);
%         end
%         saveas(h,[opt.save_folder,'combined_vol_surf_spher_spatial_profile']);
%         saveas(h,[opt.save_folder,'combined_vol_surf_spher_spatial_profile','.png']);
%     end
%     close all    
% end

end

